use crate::utils::{
    fastq_iter::open,
    position::Position,
    barcode_iter::{ validate_absolute_filepath, BarcodesIter },
    error::AppError,
    barcode_config::{ BarcodeMode, OpenstContext, validate_barcode_pattern },
};

use std::fs;
use std::io::{ Write, BufWriter };
use std::path::PathBuf;

use regex::Regex;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "filter")]
#[command(
    about = "Barcode Statistic", 
    long_about = "\
Performs comprehensive barcode distribution analysis and abundance estimation. \
This tool is mandatory when using the 'correct_bc' module, as it provides the \
essential sampling-based posterior correction counts required to accurately \
recover true barcode signals from sequencing noise."
)]
#[command(next_line_help = true)]
pub struct StatsBCArgs {
    /// The path to FASTQ file (barcode)
    #[arg(value_parser = validate_absolute_filepath)]
    input: PathBuf,

    /// Output txt file
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,

    /// barcode parsing mode
    #[arg(short = 'm', long, value_enum, default_value_t = BarcodeMode::Openst)]
    mode: BarcodeMode,

    /// Custom barcode position (only effective when mode=custom)
    ///
    /// Format: "read{1/2}:{+/-}:start-end"
    ///
    /// Due to single-ended sequencing, there should only be read1, (e.g. "read1:+:1-16" or "read1:-:2-30")
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

impl StatsBCArgs {
    pub fn init(self) -> Result<InitStatsBCArgs, AppError> {

        let output = self.output.unwrap_or_else(|| {
            let path_str = self.input.to_str().expect("Invalid UTF-8 path");
            
            let re = Regex::new(r"(?i)\.(fastq\.gz|fastq|fq\.gz|fq)$").unwrap();
            
            if re.is_match(path_str) {
                PathBuf::from(re.replace(path_str, ".txt").to_string())
            } else {
                self.input.with_extension("txt")
            }
        });

        let (pos, pattern) = match (self.mode, self.barcode_pos, self.barcode_pattern) {
            (BarcodeMode::Custom, Some(pos), Some(pattern)) => (pos, pattern),
            (BarcodeMode::Custom, _, _) => {
                return Err(
                    AppError::CommandError(
                        "Custom barcode mode requires --barcode-pos and --barcode-pattern".into()
                    )
                );
            }
            (BarcodeMode::SCMethyl, None, None) => BarcodeMode::sc_methyl(),
            (BarcodeMode::ChromiumATAC, None, None) => BarcodeMode::chromium_atac(),
            (BarcodeMode::Openst, None, None) => BarcodeMode::openst(OpenstContext::Tile),
            (BarcodeMode::ChromiumMRNA, None, None) => BarcodeMode::chromium_mrna(),
            _ => {
                unimplemented!("Other barcode modes are unimplemented!");
            }
        };

        Ok(InitStatsBCArgs {
            input: self.input,
            output: output,
            pos,
            pattern,
        })
    }
}

/// Internal struct used by add_cb pipeline
#[derive(Debug, Clone)]
pub struct InitStatsBCArgs {
    /// Input Fastq path
    input: PathBuf,
    /// Output txt path
    output: PathBuf,
    /// Barcode position specification
    pos: Position,
    /// Barcode pattern (regex)
    pattern: String,
}

impl InitStatsBCArgs {
    #[inline]
    fn input_path(&self) -> &PathBuf {
        &self.input
    }

    #[inline]
    fn output_path(&self) -> &PathBuf {
        &self.output
    }

    #[inline]
    fn pos(&self) -> &Position {
        &self.pos
    }

    #[inline]
    fn pattern(&self) -> &str {
        &self.pattern
    }

    pub fn run_stats_cb(&self) -> Result<(), AppError> {
        println!("Writing barcode list...");
        let barcode_map = BarcodesIter::new(
            open(self.input_path())?,
            self.pos(),
            self.pattern()
        ).stat_sample_barcodes(0)?;
        let mut barcode_writer = BufWriter::new(
            fs::OpenOptions::new().create(true).write(true).truncate(true).open(self.output_path())?
        );
        writeln!(barcode_writer, "Barcode\tCount")?;
        for (bc, count) in barcode_map {
            writeln!(barcode_writer, "{}\t{}", bc, count)?;
        }
        Ok(())
    }
}
