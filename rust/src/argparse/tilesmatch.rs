
use crate::argparse::touchbarcode::{validate_barcode_pattern};
use crate::utils::{
    fastqfile::{open, FastqReader},
    position::Position,
    barcode_iter::{validate_absolute_filepath, BarcodesIter},
    error::AppError,
};
use std::io;
use std::path::PathBuf;
use std::collections::HashSet;
use clap::{Parser, ValueEnum};
use rayon::prelude::*;
use rust_htslib::tbx::{self, Read};

pub fn is_valid_tile_id(value: &str) -> Result<u64, String> {
    let tile_id: u64 = value.parse()
        .map_err(|_| format!("`{}` is not valid integer", value))?;
    
    if VALID_TILE_IDS.contains(&tile_id) {
        Ok(tile_id)
    } else {
        Err(format!("tile_id {} is not in the valid range (valid range: 11101-42678)", tile_id))
    }
}

const VALID_TILE_IDS: [u64; 3744] = {
    // Array size: 4 × 2 × 6 × 78 = 3744
    let mut result = [0u64; 3744];
    let mut index = 0;
    
    let mut a = 1;
    while a <= 4 {
        let mut b = 1;
        while b <= 2 {
            let mut c = 1;
            while c <= 6 {
                let mut d = 1;
                while d <= 78 {
                    result[index] = a * 10000 + b * 1000 + c * 100 + d;
                    index += 1;
                    d += 1;
                }
                c += 1;
            }
            b += 1;
        }
        a += 1;
    }
    result
};

/// supported raw fastq.gz file or bam file
#[derive(Parser, Debug)]
#[command(name = "tilesmatch")]
#[command(about = 
    "Search for each tile that match the threshold", 
    long_about = None
)]
#[command(next_line_help = true)]
pub struct TilesMatchArgs {
    /// Generally Read1 fastq file
    #[arg(
        short = 'R', 
        long, 
        required = true,
    )]
    read: PathBuf,

    /// The path to the barcode file
    #[arg(
        short = 'I', 
        long, 
        required = true, 
        value_parser = validate_absolute_filepath,
    )]
    barcode_file: PathBuf,

    /// the tile id list to query
    #[arg(
        long, 
        value_delimiter = ' ',
        num_args = 1..,
        value_parser = is_valid_tile_id,
    )]
    tile_list: Option<Vec<u64>>,

    /// the number of barcodes used to query
    #[arg(short, long, default_value_t = 100_000_000)]
    num_barcode: usize,

    /// the threshold to filter tile
    #[arg(long, default_value_t = 0.1)]
    threshold: f32,

    /// turn on it to output tile id that passed threshold.
    #[arg(short, long)]
    quiet: bool,

    /// barcode/UMI parsing mode
    #[arg(short, long, value_enum, default_value_t = BarcodeMode::Openst)]
    mode: BarcodeMode,

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

    /// Custom barcode pattern (only effective when mode=custom)
    /// 
    /// Regex: ^[ATGCNRYMKSWHBVD]+$
    /// 
    /// there should only be the pattern before convert sequence into reverse complement sequence.
    /// (e.g. openst-barcode: VNBVNNVNNVNNVNNVNNVNNVNNVNNN, openst-seq: NNNBNNBNNBNNBNNBNNBNNBNNBVNB)
    #[arg(
        long, 
        required_if_eq("mode", "custom"), 
        value_parser = validate_barcode_pattern,
        value_name = "BARCODE_PATTERN",
    )]
    barcode_pattern: Option<String>,
}

impl TilesMatchArgs {
    pub fn init(self) -> Result<InitTilesMatchArgs, AppError> {
        let (pos, pattern) = match (self.barcode_pos, self.barcode_pattern) {
            (Some(pos), Some(pattern)) => (pos, pattern),
            (None, None) => BarcodeMode::openst(),
            _ => unreachable!("clap parse the error is impossible.")
        };
        let tile_list = if let Some(list) = self.tile_list {
            list
        } else {
            // 直接返回预生成的常量数组
            VALID_TILE_IDS.to_vec()
        };
        
        Ok(InitTilesMatchArgs::new(
            self.read, 
            self.barcode_file, 
            tile_list, 
            self.num_barcode, 
            self.threshold,
            self.quiet,
            pos,
            pattern,
        ))
    }
}

pub struct InitTilesMatchArgs {
    read: PathBuf,
    barcode_file: PathBuf,
    tile_list: Vec<u64>,
    num_barcode: usize,
    threshold: f32,
    quiet: bool,
    pos: Position,
    pattern: String,
}

impl InitTilesMatchArgs {
    #[inline]
    fn new(
        read: PathBuf,
        barcode_file: PathBuf,
        tile_list: Vec<u64>,
        num_barcode: usize,
        threshold: f32,
        quiet: bool,
        pos: Position,
        pattern: String,
    ) -> Self {
        Self { 
            read, 
            barcode_file, 
            tile_list, 
            num_barcode, 
            threshold, 
            quiet,
            pos, 
            pattern 
        }
    }

    #[inline]
    pub fn quiet(&self) -> bool { self.quiet }

    pub fn create_barcode_iter(&self) -> Result<BarcodesIter<HashSet<String>>, AppError> {
        let inner: FastqReader = open(&self.read)?;
        Ok(BarcodesIter::into_set(
            inner, 
            &self.pos, 
            &self.pattern, 
            HashSet::with_capacity(self.num_barcode)
        ))
    }

    pub fn search_tile(&self) -> Result<Vec<TileMatchReport>, AppError> {
        let barcode_list = self.create_barcode_iter()?.extract_sample_barcodes(self.num_barcode)?;
        self.tile_list.par_iter().map(
            |&tile_id| {
                let mut chip_reader = tbx::Reader::from_path(&self.barcode_file)?;
                let tid = chip_reader.tid(&tile_id.to_string())?;
                chip_reader.fetch(tid, 1000, 37100)?;

                let tile_list = chip_reader.records().map(
                    |record| {
                        let record = record?;
                        let record = unsafe { String::from_utf8_unchecked(record) };
                        let barcode = record.splitn(4, '\t').nth(3).ok_or(AppError::IoError(
                            io::Error::new(io::ErrorKind::InvalidData, "Invalid tile's barcode file format")
                        ))?;

                        Ok(barcode.to_string())
                    }
                ).collect::<Result<HashSet<String>, AppError>>()?;
                let passed_num = tile_list.intersection(&barcode_list).count();
                let percent = passed_num as f32 / tile_list.len() as f32;
                let pass_threshold = percent >= self.threshold;
                Ok(TileMatchReport::new(
                    tile_id, 
                    passed_num, 
                    tile_list.len(), 
                    percent, 
                    pass_threshold
                ))
            }
        ).collect::<Result<Vec<TileMatchReport>, AppError>>()
    }  
}

#[derive(ValueEnum, Clone, Copy, Debug)]
pub enum BarcodeMode {
    Openst,
    Custom,
}

pub type BarcodeConfig = (Position, String);
impl BarcodeMode {
    pub fn openst() -> BarcodeConfig {
        let pos = Position::new(false, false, 2, 30);
        // HDMI32-DraI: NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN
        // revcomp:     NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN
        let pattern: String = String::from("VNBVNNVNNVNNVNNVNNVNNVNNVNNN");
        (pos, pattern)
    }
}

pub struct TileMatchReport {
    tile_id: u64,
    passed_num: usize,
    total_num: usize,
    percent: f32,
    pass_threshold: bool,
}

impl TileMatchReport {
    #[inline]
    fn new(
        tile_id: u64, 
        passed_num: usize, 
        total_num: usize, 
        percent: f32, 
        pass_threshold: bool
    ) -> Self {
        Self {
            tile_id,
            passed_num,
            total_num,
            percent,
            pass_threshold,
        }
    }

    #[inline]
    pub fn tile_id(&self) -> u64 { self.tile_id }

    #[inline]
    pub fn pass_threshold(&self) -> bool { self.pass_threshold }
}

impl std::fmt::Display for TileMatchReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:<7}\t{:<12}\t{:<14}\t{:<11.5}\t{}",
            self.tile_id,
            self.total_num,
            self.passed_num,
            self.percent,
            if self.pass_threshold { 1 } else { 0 },
        )
    }
}