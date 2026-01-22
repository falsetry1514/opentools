// ------------------
// 1. std
// ------------------
use std::{ io, usize };
use std::path::PathBuf;
use std::collections::HashSet;

// ------------------
// 2. external crates
// ------------------
use clap::Parser;
use sysinfo::System;
use rayon::{ ThreadPoolBuilder, prelude::* };
use rust_htslib::tbx::{ Reader as TbxReader, Read };

// ------------------
// 3. inner crate
// ------------------
use crate::utils::{
    barcode_config::{ BarcodeMode, OpenstContext, validate_barcode_pattern },
    barcode_iter::{ validate_absolute_filepath, BarcodesIter },
    fastq_iter::open,
    position::Position,
    error::AppError,
    tile_ids::{ is_valid_tile_id, VALID_TILE_IDS },
};

/// supported raw fastq.gz file or bam file
#[derive(Parser, Debug)]
#[command(name = "tilesmatch")]
#[command(about = "Search for each tile that match the threshold", long_about = None)]
#[command(next_line_help = true)]
pub struct TilesMatchArgs {
    /// Generally Read1 fastq file
    #[arg(short = 'R', long, required = true)]
    read: PathBuf,

    /// The path to the barcode file
    #[arg(short = 'I', long, required = true, value_parser = validate_absolute_filepath)]
    barcode_file: PathBuf,

    /// the tile id list to query
    #[arg(long, value_delimiter = ' ', num_args = 1.., value_parser = is_valid_tile_id)]
    tile_list: Option<Vec<u64>>,

    /// the number of barcodes used to query; 0 for unlimited
    #[arg(short, long, default_value_t = 100_000_000)]
    num_barcode: usize,

    /// the threshold to filter tile
    #[arg(short = 'M', long, default_value_t = 0.1)]
    threshold: f32,

    /// The cores used for matching
    #[arg(short = '@', long)]
    cores: Option<usize>,

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
        value_name = "BARCODE_POS"
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
        value_name = "BARCODE_PATTERN"
    )]
    barcode_pattern: Option<String>,
}

impl TilesMatchArgs {
    pub fn init(self) -> Result<InitTilesMatchArgs, AppError> {
        let (pos, pattern) = match (self.barcode_pos, self.barcode_pattern) {
            (Some(pos), Some(pattern)) => (pos, pattern),
            (None, None) => BarcodeMode::openst(OpenstContext::Tile),
            _ => unreachable!("clap parse the error is impossible."),
        };

        let tile_list = if let Some(list) = self.tile_list {
            list
        } else {
            VALID_TILE_IDS.to_vec()
        };

        let upper_cores = {
            let mut sys = System::new();
            sys.refresh_cpu_all();
            sys.cpus().len().saturating_sub(2).max(1)
        };

        let n_core = self.cores.unwrap_or(upper_cores).min(upper_cores);

        Ok(
            InitTilesMatchArgs::new(
                self.read,
                self.barcode_file,
                tile_list,
                self.num_barcode,
                self.threshold,
                n_core,
                self.quiet,
                pos,
                pattern
            )
        )
    }
}

pub struct InitTilesMatchArgs {
    read: PathBuf,
    barcode_file: PathBuf,
    tile_list: Vec<u64>,
    num_barcode: usize,
    threshold: f32,
    cores: usize,
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
        cores: usize,
        quiet: bool,
        pos: Position,
        pattern: String
    ) -> Self {
        Self {
            read,
            barcode_file,
            tile_list,
            num_barcode,
            threshold,
            cores,
            quiet,
            pos,
            pattern,
        }
    }

    #[inline]
    pub fn quiet(&self) -> bool {
        self.quiet
    }

    #[inline]
    pub fn cores(&self) -> usize {
        self.cores
    }

    pub fn create_barcode_iter(&self) -> Result<BarcodesIter<impl std::io::Read>, AppError> {

        let inner = open(&self.read)?;
        // HashSet::with_capacity(self.num_barcode)
        Ok(
            BarcodesIter::new(
                inner,
                &self.pos,
                &self.pattern,
            )
        )
    }

    pub fn search_tile(&self) -> Result<Vec<TileMatchReport>, AppError> {
        let barcode_list = self.create_barcode_iter()?.extract_sample_barcodes(self.num_barcode)?;

        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cores())
            .build()
            .expect("Build thread pool failed");

        pool.install(|| {
            self.tile_list
                .par_iter()
                .map(|&tile_id| {
                    let mut chip_reader = TbxReader::from_path(&self.barcode_file)?;
                    let tid = chip_reader.tid(&tile_id.to_string())?;
                    chip_reader.fetch(tid, 1000, 37100)?;

                    let tile_list = chip_reader
                        .records()
                        .map(|record| {
                            let record = record?;
                            let record = unsafe { String::from_utf8_unchecked(record) };
                            let barcode = record
                                .splitn(4, '\t')
                                .nth(3)
                                .ok_or(
                                    AppError::Io(
                                        io::Error::new(
                                            io::ErrorKind::InvalidData,
                                            "Invalid tile's barcode file format"
                                        )
                                    )
                                )?;

                            Ok(barcode.to_string())
                        })
                        .collect::<Result<HashSet<String>, AppError>>()?;
                    let passed_num = tile_list.intersection(&barcode_list).count();
                    let percent = (passed_num as f32) / (tile_list.len() as f32);
                    let pass_threshold = percent >= self.threshold;
                    Ok(
                        TileMatchReport::new(
                            tile_id,
                            passed_num,
                            tile_list.len(),
                            percent,
                            pass_threshold
                        )
                    )
                })
                .collect::<Result<Vec<TileMatchReport>, AppError>>()
        })
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
    pub fn tile_id(&self) -> u64 {
        self.tile_id
    }

    #[inline]
    pub fn pass_threshold(&self) -> bool {
        self.pass_threshold
    }
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
            if self.pass_threshold {
                1
            } else {
                0
            }
        )
    }
}
