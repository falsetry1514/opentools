
use super::{
    fastqfile::{FastqReader, check_base_match, complement},
    position::Position,
    error::AppError,
};
use std::{collections::HashSet, sync::atomic::AtomicUsize};
use std::io::{Write, self};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use seq_io::fastq::Record;
use dashmap::DashSet;
use rayon::prelude::*;

pub fn validate_absolute_dirpath(s: &str) -> io::Result<PathBuf> {
    let mut path = Path::new(s).to_path_buf();
    if !path.is_dir() {
        return Err(
            io::Error::new(
                io::ErrorKind::NotADirectory,
                format!("{} is not a directory", s)
            )
        );
    }
    if path.is_relative() {
        path = path.canonicalize()?;
    }
    Ok(path)
}

pub fn validate_absolute_filepath(s: &str) -> io::Result<PathBuf> {
    let path = Path::new(s).to_path_buf();
    if !path.is_file() {
        return Err(
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{} is not a file", s)
            )
        );
    }
    Ok(path)
}

pub struct BarcodesIter<'a, W> {
    inner: FastqReader,
    pos: &'a Position,
    pattern: &'a str,
    writer: W,
}

impl<'a, W> BarcodesIter<'a, W> { 
    // Factory mathod
    pub fn new(
        inner: FastqReader, 
        pos: &'a Position, 
        pattern: &'a str, 
        writer: W,
    ) -> Self {
        Self {
            inner,
            pos,
            pattern,
            writer,
        }
    }

    // Associated method
    fn fail_quality_filter(qual: &[u8]) -> bool {
        let mut low_qual_count: u64 = 0;
        for &q in qual {
            if q < 53 { return true; }
            if q < 63 { low_qual_count += 1; }
        }
        low_qual_count > 2
    }

    fn fail_sequence_filter(seq: &[u8], pattern: &str) -> bool {
        seq.iter().zip(pattern.bytes()).any(|(&b, p)| check_base_match(b, p))
    }

    fn process_barcode(seq: &[u8], is_revcomp: bool) -> String {
        let barcode: Vec<u8> = if is_revcomp {
            seq.iter().rev().map(complement).collect()
        } else {
            seq.to_vec()
        };
        unsafe { String::from_utf8_unchecked(barcode) }
    }

    fn parse_id(id: &str) -> (&str, &str, &str, &str) {
        let mut parts = id.splitn(7, ':');
        match (parts.nth(3), parts.next(), parts.next(), parts.next()) {
            (Some(l), Some(t), Some(x), Some(y)) => (l, t, x, y),
            _ => unreachable!("Invalid fastq id occurs!"),
        }
    }
}

impl<'a, W> BarcodesIter<'a, W> 
where
    W: Write,
{
    // Factory mathod
    pub fn into_file(
        inner: FastqReader, 
        pos: &'a Position, 
        pattern: &'a str, 
        writer: W,
    ) -> Self {
        Self::new(inner, pos, pattern, writer)
    }

    // Public method
    pub fn extract_chip_barcodes(mut self) -> Result<Report, AppError> {
        let mut seen_positions = HashSet::new();
        let mut buffer = Vec::with_capacity(1000);

        let mut total_count: u64 = 0;
        let mut filter_seq_count: u64 = 0;
        let mut filter_qual_count: u64 = 0;
        let mut filter_dup_count: u64 = 0;
        for rec in self.inner.records() {
            let rec = rec?;
            total_count += 1;
            let (seq, qual) = (
                self.pos.safe_slice(&rec.seq),
                self.pos.safe_slice(&rec.qual),
            );
            let id = rec.id().expect("Invalid record id");
            let (lane, tile, x_pos, y_pos) = Self::parse_id(id);
            let pos_key = (x_pos.to_string(), y_pos.to_string());

            if Self::fail_quality_filter(qual) {
                filter_qual_count += 1;
                continue;
            }
            if Self::fail_sequence_filter(seq, self.pattern) {
                filter_seq_count += 1;
                continue; 
            }
            if !seen_positions.insert(pos_key) {
                filter_dup_count += 1;
                continue;
            }

            let barcode = Self::process_barcode(seq, self.pos.is_revcomp());
            buffer.push(format!("{}{}\t{}\t{}\t{}\n", lane, tile, x_pos, y_pos, barcode));
            if buffer.len() >= 1000 {
                self.writer.write_all(buffer.concat().as_bytes())?;
                buffer.clear();
            }
        }
        if !buffer.is_empty() {
            self.writer.write_all(buffer.concat().as_bytes())?;
        }
        self.writer.flush()?;
        
        Ok(Report::new(total_count, filter_qual_count, filter_seq_count, filter_dup_count))
    }
}

impl<'a> BarcodesIter<'a, HashSet<String>> {
    pub fn into_set(
        // tile_id: &'a str, 
        inner: FastqReader, 
        pos: &'a Position, 
        pattern: &'a str, 
        writer: HashSet<String>,
    ) -> Self {
        Self::new(inner, pos, pattern, writer)
    }

    pub fn extract_sample_barcodes(mut self, capacity: usize) -> Result<HashSet<String>, AppError> {
        let barcode_set = DashSet::new();
        let capacity_reached = AtomicBool::new(false);
        let unique_barcode_num = AtomicUsize::new(0);
        
        self.inner.records().par_bridge().try_for_each(
            |rec| -> Result<(), AppError> {
            if capacity_reached.load(Ordering::Relaxed) {
                return Ok(());
            }
            
            let rec = rec?;
            
            let (seq, qual) = (
                &rec.seq[self.pos.range()],
                &rec.qual[self.pos.range()],
            );
            
            if Self::fail_quality_filter(qual) || Self::fail_sequence_filter(seq, self.pattern) {
                return Ok(());
            }

            if capacity_reached.load(Ordering::Relaxed) {
                return Ok(());
            }
            
            let barcode = Self::process_barcode(seq, self.pos.is_revcomp());
            
            if barcode_set.insert(barcode) {
                let count = unique_barcode_num.fetch_add(1, Ordering::Relaxed) + 1;
                if count >= capacity {
                    capacity_reached.store(true, Ordering::Relaxed);
                }
            }
            Ok(())
        })?;
        Ok(barcode_set.into_iter().take(capacity).collect())
    }
}

pub struct Report {
    total_count: u64,
    filter_qual_count: u64,
    filter_seq_count: u64,
    filter_dup_count: u64,
}

impl Report {
    #[inline]
    fn new(
        total_count: u64, 
        filter_qual_count: u64, 
        filter_seq_count: u64, 
        filter_dup_count: u64
    ) -> Self {
        Self { total_count, filter_qual_count, filter_seq_count, filter_dup_count }
    }

    #[inline]
    fn filtered_count(&self) -> u64 {
        self.filter_qual_count + self.filter_seq_count + self.filter_dup_count
    }

    #[inline]
    fn passed_count(&self) -> u64 {
        self.total_count - self.filtered_count()
    }
}

impl std::fmt::Display for Report {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Total={}, Filtered={} (Qual={}, Seq={}, Dup={}), Passed={}",
            self.total_count,
            self.filtered_count(),
            self.filter_qual_count,
            self.filter_seq_count,
            self.filter_dup_count,
            self.passed_count()
        )
    }
}

