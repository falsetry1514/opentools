use crate::argparse::touchbarcode::TouchBarcodeReport;
use crate::utils::{
    error::AppError,
    position::Position,
};
use seq_io::fastq::{ Record, Reader};
use std::collections::HashSet;
use std::io::{ self, Write, Read };
use std::path::{ Path, PathBuf };

pub fn validate_absolute_dirpath(s: &str) -> io::Result<PathBuf> {
    let mut path = Path::new(s).to_path_buf();
    if !path.is_dir() {
        return Err(
            io::Error::new(io::ErrorKind::NotADirectory, format!("{} is not a directory", s))
        );
    }
    if path.is_relative() {
        path = path.canonicalize()?;
    }
    Ok(path)
}

pub fn validate_absolute_filepath(s: &str) -> io::Result<PathBuf> {
    let path = Path::new(s).to_path_buf();
    if path.exists() && !path.is_file() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, format!("{} is not a file", s)));
    }
    Ok(path)
}

fn complement(b: &u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => unreachable!("Invalid base: {b}"),
    }
}

pub fn check_base_match(base: u8, pattern_char: u8) -> bool {
    // 碱基匹配
    match (base, pattern_char) {
        (b'N', _) => true,
        (b'A', b'A' | b'R' | b'M' | b'W' | b'H' | b'V' | b'D' | b'N') => false,
        (b'T', b'T' | b'Y' | b'K' | b'W' | b'H' | b'B' | b'D' | b'N') => false,
        (b'G', b'G' | b'R' | b'K' | b'S' | b'B' | b'V' | b'D' | b'N') => false,
        (b'C', b'C' | b'Y' | b'M' | b'S' | b'H' | b'B' | b'V' | b'N') => false,
        _ => true,
    }
}

fn parse_id(id: &str) -> (&str, &str, &str, &str) {
    let mut parts = id.splitn(7, ':');
    match (parts.nth(3), parts.next(), parts.next(), parts.next()) {
        (Some(l), Some(t), Some(x), Some(y)) => (l, t, x, y),
        _ => unreachable!("Invalid fastq id occurs!"),
    }
}

pub fn process_sequence(seq: &[u8], revcomp: bool) -> String {
    let seq = if revcomp {
        seq.iter().rev().map(complement).collect()
    } else {
        seq.to_vec()
    };
    unsafe { String::from_utf8_unchecked(seq) }
}

pub struct BarcodesIter<'a, R: Read> {
    inner: Reader<R>,
    pos: &'a Position,
    pattern: &'a str,
}

impl<'a, R: Read> BarcodesIter<'a, R> {
    // Factory mathod
    #[inline]
    pub fn new(inner: Reader<R>, pos: &'a Position, pattern: &'a str) -> Self {
        Self {
            inner,
            pos,
            pattern,
        }
    }

    #[inline]
    pub fn pos(&self) -> &Position {
        &self.pos
    }

    pub fn inner(&mut self) -> &mut Reader<R> {
        &mut self.inner
    }

    fn for_each_barcode<F>(
        inner: &mut Reader<R>,
        pos: &Position,
        mut f: F
    ) -> Result<(), AppError>
        where F: FnMut(BarcodeRecord) -> Result<(), AppError>
    {
        for rec in inner.records() {
            let rec = rec?;
            let id = rec.id().expect("Invalid record id");
            let seq = pos.safe_slice(&rec.seq);
            let qual = pos.safe_slice(&rec.qual);
            let (lane, tile, x_pos, y_pos) = parse_id(id);
            f(BarcodeRecord::new(seq, qual, lane, tile, x_pos, y_pos, pos.is_revcomp()))?;
        }
        Ok(())
    }

    // Public method
    pub fn extract_chip_barcodes<W: Write>(mut self, mut writer: W) -> Result<TouchBarcodeReport, AppError> {
        let mut seen_positions = HashSet::new();
        let mut buffer = Vec::with_capacity(1000);

        let mut total_count: u64 = 0;
        let mut filter_seq_count: u64 = 0;
        let mut filter_qual_count: u64 = 0;
        let mut filter_dup_count: u64 = 0;

        Self::for_each_barcode(&mut self.inner, self.pos, |br| {
            total_count += 1;
            let pos_key = (br.x_pos().to_string(), br.y_pos().to_string());

            if br.fail_quality_filter() {
                filter_qual_count += 1;
                return Ok(());
            }
            if br.fail_sequence_filter(self.pattern) {
                filter_seq_count += 1;
                return Ok(());
            }
            if !seen_positions.insert(pos_key) {
                filter_dup_count += 1;
                return Ok(());
            }

            buffer.push(
                format!("{}{}\t{}\t{}\t{}\n", br.lane(), br.tile(), br.x_pos(), br.y_pos(), br.process_sequence())
            );

            if buffer.len() >= 1000 {
                for line in &buffer {
                    writer.write_all(line.as_bytes())?;
                }
                buffer.clear();
            }

            Ok(())
        })?;

        if !buffer.is_empty() {
            for line in &buffer {
                writer.write_all(line.as_bytes())?;
            }
        }
        writer.flush()?;

        Ok(TouchBarcodeReport::new(total_count, filter_qual_count, filter_seq_count, filter_dup_count))
    }

    pub fn extract_sample_barcodes(mut self, capacity: usize) -> Result<HashSet<String>, AppError> {
        let mut barcode_set = HashSet::new();
        let mut unique_barcode_num = 0;

        Self::for_each_barcode(&mut self.inner, self.pos, |br| {
            if barcode_set.insert(br.process_sequence()) {
                unique_barcode_num += 1;
                if capacity != 0 && unique_barcode_num >= capacity {
                    // 直接提前终止
                    return Err(AppError::EarlyStop);
                }
            }
            Ok(())
        }).or_else(|e| {
            if matches!(e, AppError::EarlyStop) { Ok(()) } else { Err(e) }
        })?;

        Ok(barcode_set)
    }
}

impl<'a, R: Read> Iterator for BarcodesIter<'a, R> {
    type Item = Result<String, AppError>;

    fn next(&mut self) -> Option<Self::Item> {
        let rec = match self.inner.next() {
            Some(Ok(bc)) => bc.to_owned_record(),
            Some(Err(e)) => return Some(Err(e.into())),
            None => return None,
        };

        let seq = self.pos.safe_slice(&rec.seq);
        let qual = self.pos.safe_slice(&rec.qual);
        let id = rec.id().expect("Invalid record id");
        let (lane, tile, x_pos, y_pos) = parse_id(id);

        let barcode = BarcodeRecord::new(seq, qual, lane, tile, x_pos, y_pos, self.pos.is_revcomp()).process_sequence();

        Some(Ok(barcode))
    }
}

pub struct BarcodeRecord<'a> {
    seq: &'a [u8],
    qual: &'a [u8],
    lane: &'a str,
    tile: &'a str,
    x: &'a str,
    y: &'a str,
    revcomp: bool,
}

impl<'a> BarcodeRecord<'a> {
    #[inline]
    pub fn new(
        seq: &'a [u8],
        qual: &'a [u8],
        lane: &'a str,
        tile: &'a str,
        x: &'a str,
        y: &'a str,
        revcomp: bool
    ) -> Self {
        Self {
            seq,
            qual,
            lane,
            tile,
            x,
            y,
            revcomp,
        }
    }

    #[inline]
    pub fn seq(&self) -> &[u8] {
        self.seq
    }

    #[inline]
    pub fn qual(&self) -> &[u8] {
        self.qual
    }

    #[inline]
    pub fn lane(&self) -> &str {
        self.lane
    }

    #[inline]
    pub fn tile(&self) -> &str {
        self.tile
    }

    #[inline]
    pub fn x_pos(&self) -> &str {
        self.x
    }

    #[inline]
    pub fn y_pos(&self) -> &str {
        self.y
    }

    #[inline]
    pub fn revcomp(&self) -> bool {
        self.revcomp
    }

    // Associated method
    pub fn fail_quality_filter(&self) -> bool {
        let mut low_qual_count: u64 = 0;
        for &q in self.qual {
            if q < 53 {
                return true;
            }
            if q < 63 {
                low_qual_count += 1;
            }
        }
        low_qual_count > 2
    }

    pub fn fail_sequence_filter(&self, pattern: &str) -> bool {
        self.seq
            .iter()
            .zip(pattern.bytes())
            .any(|(&b, p)| check_base_match(b, p))
    }

    pub fn process_sequence(&self) -> String {
        let seq = if self.revcomp {
            self.seq.iter().rev().map(complement).collect()
        } else {
            self.seq.to_vec()
        };
        unsafe { String::from_utf8_unchecked(seq) }
    }
}


