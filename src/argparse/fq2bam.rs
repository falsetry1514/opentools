
use crate::utils::{
    position::Position, 
    fastqfile::{FastqReader, open, replace_asterisk, complement},
    error::AppError,
};

use std::io;
use std::path::PathBuf;
use clap::{Parser, ValueEnum};
use seq_io::fastq::{Record, RefRecord};
use rust_htslib::{bam, bam::header::HeaderRecord, bam::Record as BamRecord, bam::record::Aux};

const VERSION: &str = &env!("CARGO_PKG_VERSION");

// fq2bam subcommand
#[derive(Parser, Debug)]
#[command(name = "fq2bam")]
#[command(about = 
    "Convert paired PE150 fastq.gz into bam file in a specfied mode", 
    long_about = None
)]
#[command(next_line_help = true)]
pub struct Fq2BamArgs {
    ///Read1 fastq files (comma separated)
    #[arg(
        short = '1', 
        long, 
        required = true, 
        value_delimiter = ',',
    )]
    read1: Vec<PathBuf>,

    /// Read2 fastq files (comma separated)
    #[arg(
        short = '2', 
        long, 
        required = true, 
        value_delimiter = ',',
    )]
    read2: Vec<PathBuf>,

    /// sample info in bam header [default: Unknown]
    #[arg(long, default_value = "Unknown")]
    sample_id: String,

    /// output format [default: bam]
    #[arg(long, value_enum, default_value_t = Format::Bam)]
    format: Format,

    /// barcode/UMI parsing mode
    #[arg(short, long, value_enum, default_value_t = Mode::Openst)]
    mode: Mode,

    /// Custom barcode position (only effective when mode=custom)
    /// 
    /// Format: "read{1/2}:{+/-}:start-end" 
    /// 
    /// (e.g. "read1:+:1-16" or "read2:-:20-end")
    #[arg(
        long, 
        required_if_eq("mode", "custom"), 
        value_parser = clap::value_parser!(Position),
        value_name = "BARCODE_POS",
    )]
    barcode_pos: Option<Position>,

    /// Custom UMI position (only effective when mode=custom)
    /// 
    /// Format: "read{1/2}:{+/-}:start-end" 
    /// 
    /// (e.g. "read1:+:1-16" or "read2:-:20-end")
    #[arg(
        long, 
        required_if_eq("mode", "custom"), 
        value_parser = clap::value_parser!(Position),
        value_name = "UMI_POS",
    )]
    umi_pos: Option<Position>,

    /// Custom Read position (only effective when mode=custom)
    /// 
    /// Format: "read{1/2}:{+/-}:start-end" 
    /// 
    /// (e.g. "read1:+:1-16" or "read2:-:20-end")
    #[arg(
        long, 
        required_if_eq("mode", "custom"), 
        value_parser = clap::value_parser!(Position),
        value_name = "READ_POS",
    )]
    read_pos: Option<Position>,
}

impl Fq2BamArgs {

    #[inline]
    pub fn format(&self) -> &Format { &self.format }

    #[inline]
    pub fn read_pos(&self) -> Option<&Position> { self.read_pos.as_ref() }

    #[inline]
    pub fn barcode_pos(&self) -> Option<&Position> { self.barcode_pos.as_ref() }

    #[inline]
    pub fn umi_pos(&self) -> Option<&Position> { self.umi_pos.as_ref() }

    #[inline]
    pub fn validate_eq_file_count(&self) -> Result<(), AppError> {
        if self.read1.len() != self.read2.len() {
            return Err(
                AppError::NotEqualFileNumber { 
                    n1: self.read1.len(), n2: self.read2.len()
                }
            )
        }
        Ok(())
    }

    pub fn create_bam_header(&self) -> bam::Writer {
        let mut header = bam::Header::new();
        header.push_record(&HeaderRecord::new(b"HD\tVN:1.6\tSO:unsorted"));
        header.push_record(&HeaderRecord::new(
            format!("RG\tID:A\tSM:{}", self.sample_id).as_bytes()
        ));
        header.push_record(&HeaderRecord::new(
            format!(
                "PG\tPN:opentools\tID:fq2bam\tVN:{}\tCL:{}",
                VERSION,
                std::env::args().collect::<Vec<String>>().join(" "),
            ).as_bytes(),
        ));

        match self.format {
            Format::Bam => bam::Writer::from_stdout(&mut header, bam::Format::Bam).unwrap(),
            Format::Sam => bam::Writer::from_stdout(&mut header, bam::Format::Sam).unwrap(),
        }
    }

    pub fn paired_readers(&self) -> impl Iterator<Item = io::Result<PairedFastqReader>> {
        self.read1.iter()
            .zip(self.read2.iter())
            .map(
                |(r1, r2)| 
                PairedFastqReader::from_paths(r1, r2)
            )
    }

    pub fn record_config(&mut self) -> BamConfig {
        match self.mode {
            Mode::Openst => Mode::openst(),
            Mode::OpenTSO => Mode::open_tso(),
            Mode::Custom => {
                Mode::custom(
                    self.barcode_pos.take().unwrap(), 
                    self.umi_pos.take().unwrap(), 
                    self.read_pos.take().unwrap(),
                )
            },
        }
    }
}



#[derive(ValueEnum, Clone, Copy, Debug)]
pub enum Format {
    Bam,
    Sam,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum Mode {
    Openst, // 使用预定义位置
    OpenTSO,
    Custom,
}

impl Mode {

    #[inline]
    fn openst() -> BamConfig {
        BamConfig::new(
            Position::new(false, false, 2, 30), 
            Position::new(true, false, 0, 9), 
            Position::new(true, false, 9, 150),
        )
    }

    #[inline]
    fn open_tso() -> BamConfig {
        BamConfig::new(
            Position::new(false, false, 2, 30), 
            Position::new(false, false, 12, 20), 
            Position::new(true, false, 9, 150),
        )
    }

    #[inline]
    fn custom(cr: Position, ur: Position, read: Position) -> BamConfig {
        BamConfig::new(cr, ur, read)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct BamConfig {
    barcode_pos: Position,
    umi_pos: Position,
    read_pos: Position,
}

impl BamConfig {
    pub fn new(barcode_pos: Position, umi_pos: Position, read_pos: Position) -> Self {
        Self {
            barcode_pos, umi_pos, read_pos,
        }
    }

    #[inline]
    pub fn barcode_pos(&self) -> &Position { &self.barcode_pos }

    #[inline]
    pub fn umi_pos(&self) -> &Position { &self.umi_pos }

    #[inline]
    pub fn read_pos(&self) -> &Position { &self.read_pos }
}

impl std::fmt::Display for BamConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f, "\n\tBarcode: {}\n\tUMI: {}\n\tRead: {}", 
            self.barcode_pos, self.umi_pos, self.read_pos,
        )
    }
}

pub struct PairedFastqReader {
    reader1: FastqReader,
    reader2: FastqReader,
}

impl PairedFastqReader {

    #[inline]
    pub fn from_paths(p1: &PathBuf, p2: &PathBuf) -> io::Result<Self> {
        let reader1: FastqReader = open(p1)?;
        let reader2: FastqReader = open(p2)?;
        Ok(Self { reader1, reader2 })
    }

    #[inline]
    pub fn records(&mut self, config: BamConfig) -> RecordsIter {
        RecordsIter { inner: self, config }
    }

    fn next_record(&mut self) -> Option<Result<PairedOwnedRecord, AppError>> {
        let r1 = self.reader1.next();
        let r2 = self.reader2.next();

        match (r1, r2) {
            (None, None) => None,
            (None, _) | (_, None) => {
                Some(Err(
                    AppError::FastqFileLengthNotEqual {
                        line1: self.reader1.position().line(),
                        line2: self.reader2.position().line(),
                    }
                ))
            },
            (Some(Ok(rec1)), Some(Ok(rec2))) 
            if rec1.id_bytes() == rec2.id_bytes() => {
                Some(Ok(PairedOwnedRecord::new(&rec1, &rec2)))
            },
            (Some(Ok(rec1)), Some(Ok(rec2))) => {
                let id1 = rec1.id().unwrap_or_default().to_owned();
                let id2 = rec2.id().unwrap_or_default().to_owned();
                let line = self.reader1.position().line();
                Some(Err(
                    AppError::FastqIdMismatch { line, id1, id2 }
                ))
            },
            (Some(Err(e)), Some(_)) | (Some(Ok(_)), Some(Err(e))) => {
                Some(Err(AppError::FastqParseError(e)))
            },
        }
    }
}

pub struct RecordsIter<'a> {
    inner: &'a mut PairedFastqReader,
    config: BamConfig,
}

impl<'a> Iterator for RecordsIter<'a> {
    type Item = Result<BamRecord, AppError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next_record() {
            None => None,
            Some(Err(e)) => Some(Err(e)),
            Some(Ok(record)) => match self.parse_and_validate_record(record) {
                Ok(bam_record) => Some(Ok(bam_record)),
                Err(e) => Some(Err(e)),
            }
        }
    }
}

impl<'a> RecordsIter<'a> {
    fn parse_and_validate_record(&self, record: PairedOwnedRecord) -> Result<BamRecord, AppError> {        
        // 解析主读取区域
        let (seq, qual) = self.parse_read_area(&record);
        // 解析条形码和UMI
        let (cr, cy) = self.parse_tag_area(&record, &self.config.barcode_pos());
        let (ur, uy) = self.parse_tag_area(&record, &self.config.umi_pos());

        // 构建BAM记录
        let mut bam_record = BamRecord::new();
        bam_record.set(&record.qname, None, &seq, &qual);
        bam_record.push_aux(b"CR", Aux::String(&String::from_utf8_lossy(&cr)))?;
        bam_record.push_aux(b"CY", Aux::String(&String::from_utf8_lossy(&cy)))?;
        bam_record.push_aux(b"UR", Aux::String(&String::from_utf8_lossy(&ur)))?;
        bam_record.push_aux(b"UY", Aux::String(&String::from_utf8_lossy(&uy)))?;
        
        Ok(bam_record)
    }

    fn parse_read_area(&self, record: &PairedOwnedRecord) -> (Vec<u8>, Vec<u8>) {
        let pos = self.config.read_pos();
        let (seq, qual) = if pos.is_read2() {
            (pos.safe_slice(&record.r2), pos.safe_slice(&record.q2))
        } else {
            (pos.safe_slice(&record.r1), pos.safe_slice(&record.q1))
        };
        
        if pos.is_revcomp() {
            (
                seq.iter().rev().map(|b| complement(b)).collect(),
                qual.iter().rev().map(|q| replace_asterisk(q) - 33).collect()
            )
        } else {
            (
                seq.to_vec(),
                qual.iter().map(|q| replace_asterisk(q) - 33).collect()
            )
        }
    }

    fn parse_tag_area(&self, record: &PairedOwnedRecord, pos: &Position) -> (Vec<u8>, Vec<u8>) {
        let (seq_data, qual_data) = if pos.is_read2() {
            (&record.r2[pos.range()], &record.q2[pos.range()])
        } else {
            (&record.r1[pos.range()], &record.q1[pos.range()])
        };
        if pos.is_revcomp() {
            (
                seq_data.iter().rev().map(|b| complement(b)).collect(),
                qual_data.iter().rev().map(|q| replace_asterisk(q)).collect(),
            )
        } else {
            (
                seq_data.to_vec(),
                qual_data.iter().map(|q| replace_asterisk(q)).collect()
            )
        }
    }
}

#[derive(Debug)]
pub struct PairedOwnedRecord {
    qname: Vec<u8>,
    r1: Vec<u8>,
    r2: Vec<u8>,
    q1: Vec<u8>,
    q2: Vec<u8>,
}

impl PairedOwnedRecord {
    #[inline]
    pub fn new(r1: &RefRecord, r2: &RefRecord) -> Self {
        Self {
            qname: r1.id_bytes().to_vec(),
            r1: r1.seq().to_vec(),
            r2: r2.seq().to_vec(),
            q1: r1.qual().to_vec(),
            q2: r2.qual().to_vec(),
        }
    }
}