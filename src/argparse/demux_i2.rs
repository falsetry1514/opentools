use crate::utils::{
    fastq_iter::open,
    position::Position,
    barcode_iter::{
        validate_absolute_filepath,
        validate_absolute_dirpath,
        BarcodesIter,
        process_sequence,
    },
    error::AppError,
    barcode_config::{ BarcodeMode, validate_barcode_pattern },
    demux_iter::{ DemuxIter, DemuxWriters },
};

use std::{ fs::{ self, File }, io, path::Path };
use std::io::{ Write, BufWriter };
use std::path::PathBuf;
use std::collections::{ HashMap, HashSet };

use sysinfo::System;

use flate2::Compression;
use flate2::write::GzEncoder;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "demux-i2")]
#[command(about = "Demultiplex R1/R2 FASTQ by I2 barcode sequence", long_about = None)]
#[command(next_line_help = true)]
pub struct DemuxI2Args {
    /// The path to Read 1 FASTQ file (gz supported)
    #[arg(short = '1', long, required = true, value_parser = validate_absolute_filepath)]
    r1: PathBuf,

    /// The path to Read 2 FASTQ file (gz supported)
    #[arg(short = '2', long, required = true, value_parser = validate_absolute_filepath)]
    r2: PathBuf,

    /// The path to Index 2 FASTQ file (barcode)
    #[arg(short = 'I', long, required = true, value_parser = validate_absolute_filepath)]
    i2: PathBuf,

    /// Output directory for demultiplexed FASTQ files
    #[arg(
        short,
        long,
        required = true,
        value_parser = validate_absolute_dirpath,
        default_value_os = "."
    )]
    output_dir: PathBuf,

    /// Temporary directory used for bucketed processing
    #[arg(long, value_parser = validate_absolute_dirpath)]
    tmp_dir: Option<PathBuf>,

    /// Number of threads used for demultiplexing
    #[arg(short = '@', long)]
    cores: Option<usize>,

    /// Barcode length used for Stage 1 bucketing, n(total bucket) = 3^(bucket_len)
    #[arg(short = 'p', long, default_value_t = 6)]
    bucket_len: usize,

    /// Write barcode into read name and CB tag
    #[arg(long, default_value_t = true)]
    emit_cb: bool,

    /// barcode parsing mode
    #[arg(short, long, value_enum, default_value_t = BarcodeMode::SCMethyl)]
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

impl DemuxI2Args {
    pub fn init(self) -> Result<InitDemuxI2Args, AppError> {
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
            _ => {
                unimplemented!("Other barcode modes are unimplemented!");
            }
        };

        let tmp_dir = self.tmp_dir.unwrap_or(self.output_dir.join("tmp"));

        let upper_cores = {
            let mut sys = System::new();
            sys.refresh_cpu_all();
            sys.cpus().len().saturating_sub(2).max(1)
        };

        let n_core = self.cores.unwrap_or(upper_cores).min(upper_cores);

        Ok(InitDemuxI2Args {
            r1: self.r1,
            r2: self.r2,
            i2: self.i2,
            output_dir: self.output_dir,
            tmp_dir,
            cores: n_core,
            bucket_len: self.bucket_len,
            emit_cb: self.emit_cb,
            pos,
            pattern,
        })
    }
}

/// Internal struct used by demux pipeline
/// Holds all resolved parameters ready for Stage1/Stage2 processing
#[derive(Debug, Clone)]
pub struct InitDemuxI2Args {
    /// Path to Read1 FASTQ
    r1: PathBuf,
    /// Path to Read2 FASTQ
    r2: PathBuf,
    /// Path to Index2 FASTQ
    i2: PathBuf,
    /// Output directory for demultiplexed FASTQ
    output_dir: PathBuf,
    /// Temporary directory for Stage1 bucketing
    tmp_dir: PathBuf,
    /// Number of threads to use (capped by system cores)
    cores: usize,
    /// Number of bases from barcode prefix for Stage1 bucketing
    bucket_len: usize,
    /// Whether to write barcode into read name / CB tag
    emit_cb: bool,
    /// Barcode position specification
    pos: Position,
    /// Barcode pattern (regex)
    pattern: String,
}

impl InitDemuxI2Args {
    #[inline]
    pub fn r1(&self) -> &Path {
        self.r1.as_path()
    }

    #[inline]
    pub fn r2(&self) -> &Path {
        self.r2.as_path()
    }

    #[inline]
    pub fn i2(&self) -> &Path {
        self.i2.as_path()
    }

    #[inline]
    pub fn output_dir(&self) -> &Path {
        self.output_dir.as_path()
    }

    #[inline]
    pub fn tmp_dir(&self) -> &Path {
        self.tmp_dir.as_path()
    }

    #[inline]
    pub fn cores(&self) -> usize {
        self.cores
    }

    #[inline]
    pub fn bucket_len(&self) -> usize {
        self.bucket_len
    }

    #[inline]
    pub fn emit_cb(&self) -> bool {
        self.emit_cb
    }

    #[inline]
    fn pos(&self) -> &Position {
        &self.pos
    }

    #[inline]
    fn pattern(&self) -> &str {
        &self.pattern
    }

    pub fn get_barcode(&self, seq: &[u8]) -> String {
        process_sequence(self.pos.safe_slice(seq), self.pos().is_revcomp())
    }

    pub fn create_demuxread_iter<P>(&self, r1: P, r2: P, i2: P) -> io::Result<DemuxIter>
        where P: AsRef<Path>
    {
        Ok(
            DemuxIter::new(
                open(r1)?,
                open(r2)?,
                BarcodesIter::new(open(i2)?, self.pos(), self.pattern())
            )
        )
    }

    pub fn render_writer(&self, bucket: &str) -> io::Result<DemuxWriters> {
        let bucket_dir = self.tmp_dir().join(bucket);
        fs::create_dir_all(&bucket_dir)?;

        let r1 = BufWriter::new(
            fs::OpenOptions::new().write(true).create(true).open(bucket_dir.join("R1.fastq"))?
        );

        let r2 = BufWriter::new(
            fs::OpenOptions::new().write(true).create(true).open(bucket_dir.join("R2.fastq"))?
        );

        let i2 = BufWriter::new(
            fs::OpenOptions::new().write(true).create(true).open(bucket_dir.join("I2.fastq"))?
        );

        Ok(DemuxWriters { r1, r2, i2 })
    }

    pub fn stage1_bucket(&self) -> Result<HashSet<String>, AppError> {
        let mut bucket_map: HashMap<String, DemuxWriters> = HashMap::new();

        for demux_rec in self.create_demuxread_iter(self.r1(), self.r2(), self.i2())? {
            let demux_rec = demux_rec?;

            let bc = self.get_barcode(demux_rec.barcode().as_bytes());
            let bucket = bc.chars().take(self.bucket_len()).collect::<String>();

            let writers = bucket_map
                .entry(bucket.clone())
                .or_insert_with(|| {
                    self.render_writer(&bucket).expect("failed to create bucket writer")
                });

            demux_rec.write_r1(&mut writers.r1)?;
            demux_rec.write_r2(&mut writers.r2)?;
            demux_rec.write_i2(&mut writers.i2)?;
        }

        Ok(bucket_map.into_keys().collect::<HashSet<String>>())
    }

    pub fn stage2_split_fastq(
        &self,
        buckets: HashSet<String>
    ) -> Result<(), AppError> {
        let out_dir = self.output_dir().join("split_fastq");
        fs::create_dir_all(&out_dir)?;

        for bucket in buckets {
            let bucket_path = self.tmp_dir().join(&bucket);
            if !bucket_path.is_dir() {
                continue;
            }

            let mut writer_map: HashMap<
                String,
                (BufWriter<GzEncoder<File>>, BufWriter<GzEncoder<File>>)
            > = HashMap::new();

            for demux_rec in self.create_demuxread_iter(
                bucket_path.join("R1.fastq").as_path(),
                bucket_path.join("R2.fastq").as_path(),
                bucket_path.join("I2.fastq").as_path()
            )? {
                let demux_rec = demux_rec?;

                let bc = self.get_barcode(demux_rec.barcode().as_bytes());
                let bc_dir = out_dir.join(&bc);

                let writers = writer_map.entry(bc).or_insert_with(|| {
                    fs::create_dir_all(&bc_dir).unwrap();

                    let r1_file = File::create(bc_dir.join("cell_R1.fastq.gz")).unwrap();
                    let r2_file = File::create(bc_dir.join("cell_R2.fastq.gz")).unwrap();

                    (
                        BufWriter::new(GzEncoder::new(r1_file, Compression::default())),
                        BufWriter::new(GzEncoder::new(r2_file, Compression::default())),
                    )
                });

                demux_rec.write_r1(&mut writers.0)?;
                demux_rec.write_r2(&mut writers.1)?;
            }

            for (mut r1, mut r2) in writer_map.into_values() {
                r1.flush()?;
                r2.flush()?;
            }
        }
        Ok(())
    }

    pub fn run_demux(&self) -> Result<(), AppError> {
        if !self.output_dir().exists() {
            fs::create_dir_all(self.output_dir())?;
        }
        fs::create_dir_all(&self.tmp_dir())?;

        println!("Stage 1: Bucketing reads by barcode prefix...");
        let buckets = self.stage1_bucket()?; // Stage1仍可串行或并行

        println!("Stage 2: Demultiplexing within buckets...");
        self.stage2_split_fastq(buckets)?; // Stage2并行

        println!("Cleaning up temporary files...");
        fs::remove_dir_all(&self.tmp_dir)?;

        Ok(())
    }
}
