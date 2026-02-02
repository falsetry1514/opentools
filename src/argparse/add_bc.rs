use crate::utils::{
    fastq_iter::open,
    position::Position,
    barcode_iter::{ validate_absolute_filepath, BarcodesIter },
    error::AppError,
    barcode_config::{ BarcodeMode, OpenstContext, validate_barcode_pattern },
};

use std::io;
use std::path::PathBuf;

use rust_htslib::bam::{ self, Read };
use rust_htslib::bam::record::Aux;
use clap::Parser;
use itertools::Itertools;

#[derive(Parser, Debug)]
#[command(name = "addbc")]
#[command(about = "Add I2 barcode sequence to BAM tag", long_about = None)]
#[command(next_line_help = true)]
pub struct AddBCArgs {
    /// Input BAM file (use '-' for stdin)
    #[arg()]
    input: String,

    /// Output BAM file (use '-' for stdout)
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// The path to Index 2 FASTQ file (barcode)
    #[arg(short = 'I', long, required = true, value_parser = validate_absolute_filepath)]
    i2: PathBuf,

    /// barcode parsing mode
    #[arg(short = 'm', long, value_enum, default_value_t = BarcodeMode::SCMethyl)]
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

impl AddBCArgs {
    pub fn init(self) -> Result<InitAddBCArgs, AppError> {
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

        Ok(InitAddBCArgs {
            input: self.input,
            output: self.output,
            i2: self.i2,
            pos,
            pattern,
        })
    }
}

/// Internal struct used by add_cb pipeline
#[derive(Debug, Clone)]
pub struct InitAddBCArgs {
    /// Input BAM path (can be "-" for stdin)
    input: String,
    /// Output BAM path (can be "-" for stdout)
    output: Option<String>,
    /// Path to Index2 FASTQ
    i2: PathBuf,
    /// Barcode position specification
    pos: Position,
    /// Barcode pattern (regex)
    pattern: String,
}

impl InitAddBCArgs {
    #[inline]
    fn pos(&self) -> &Position {
        &self.pos
    }

    #[inline]
    fn pattern(&self) -> &str {
        &self.pattern
    }

    /// Main function to add barcode to BAM - SIMPLIFIED VERSION
    /// Assumes BAM records are in same order as original FASTQ (interleaved: read1, read2, read1, read2, ...)
    pub fn run_add_cb(&self) -> Result<(), AppError> {
        eprintln!("Starting to add barcodes to BAM...");
        eprintln!("Assumption: BAM file preserves original FASTQ order (interleaved)");

        // open input bam
        let mut in_bam = if self.input == "-" {
            bam::Reader
                ::from_stdin()
                .map_err(|e| {
                    AppError::Io(
                        io::Error::new(
                            io::ErrorKind::Other,
                            format!("Failed to read BAM from stdin: {}", e)
                        )
                    )
                })?
        } else {
            bam::Reader
                ::from_path(&self.input)
                .map_err(|e| {
                    AppError::Io(
                        io::Error::new(
                            io::ErrorKind::Other,
                            format!("Failed to open BAM file {}: {}", self.input, e)
                        )
                    )
                })?
        };

        // 获取header并创建输出BAM
        let header = bam::Header::from_template(in_bam.header());
        let mut out_bam = if let Some(output) = &self.output {
            if output == "-" {
                bam::Writer
                    ::from_stdout(&header, bam::Format::Bam)
                    .map_err(|e| {
                        AppError::Io(
                            io::Error::new(
                                io::ErrorKind::Other,
                                format!("Failed to write BAM to stdout: {}", e)
                            )
                        )
                    })?
            } else {
                bam::Writer
                    ::from_path(output, &header, bam::Format::Bam)
                    .map_err(|e| {
                        AppError::Io(
                            io::Error::new(
                                io::ErrorKind::Other,
                                format!("Failed to create output BAM file {}: {}", output, e)
                            )
                        )
                    })?
            }
        } else {
            bam::Writer
                ::from_stdout(&header, bam::Format::Bam)
                .map_err(|e| {
                    AppError::Io(
                        io::Error::new(
                            io::ErrorKind::Other,
                            format!("Failed to write BAM to stdout: {}", e)
                        )
                    )
                })?
        };

        let i2_fastq = open(&self.i2)?;
        let barcode_iter = BarcodesIter::new(i2_fastq, self.pos(), self.pattern());

        let primary_bam_iter = in_bam.records().filter_map(|rec_res| {
            match rec_res {
                Ok(rec) if !rec.is_secondary() && !rec.is_supplementary() => Some(rec),
                _ => None,
            }
        });

        let mut primary_bam_iter = primary_bam_iter.peekable();
        let mut barcode_iter = barcode_iter.peekable();

        for ((mut record1, mut record2), barcode) in primary_bam_iter
            .by_ref()
            .tuples()
            .zip(barcode_iter.by_ref()) {
            let barcode = barcode?;

            record1.push_aux(b"CR", Aux::String(&barcode))?;
            record2.push_aux(b"CR", Aux::String(&barcode))?;

            out_bam.write(&record1)?;
            out_bam.write(&record2)?;
        }

        if primary_bam_iter.peek().is_some() {
            eprintln!("Opentools Add Barcode - Warning ⚠️ - BAM has extra reads not matched by barcodes");
        }
        if barcode_iter.peek().is_some() {
            eprintln!("Opentools Add Barcode - Warning ⚠️ - More barcodes than BAM fragments");
        }

        Ok(())
    }
}
