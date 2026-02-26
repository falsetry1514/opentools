// ------------------
// 1. std
// ------------------
use std::fs;
use std::io::{ self, Write, BufWriter };
use std::path::PathBuf;

// ------------------
// 2. external crates
// ------------------
use clap::Parser;
use sysinfo::System;
use dashmap::DashSet;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use indicatif::{ ProgressBar, ProgressStyle };
use rust_htslib::tbx::{ Reader as TbxReader, Read };

// ------------------
// 3. inner crate
// ------------------
use crate::utils::{
    barcode_iter::validate_absolute_filepath,
    error::AppError,
    tile_ids::is_valid_tile_id,
};

#[derive(Parser, Debug)]
#[command(name = "dedupbarcode")]
#[command(about = "Deduplicate barcodes among the given tiles id list", long_about = None)]
#[command(next_line_help = true)]
pub struct DedupBarcodeArgs {
    /// The path to the barcode file
    #[arg(value_parser = validate_absolute_filepath)]
    barcode_file: PathBuf,

    /// the tile id list to query
    #[arg(long, value_delimiter = ' ', num_args = 1.., value_parser = is_valid_tile_id)]
    tile_list: Vec<u64>,

    /// The path to the FASTQ file
    #[arg(short, long, required = true)]
    output_dir: PathBuf,

    /// The cores used for deduplate
    #[arg(short = '@', long)]
    cores: Option<usize>,
}

impl DedupBarcodeArgs {
    pub fn init(self) -> Result<InitDedupBarcodeArgs, AppError> {
        let parent = self.output_dir.parent().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "output_dir has no parent directory",
            )
        })?;
        let canon_parent = match parent.canonicalize() {
            Ok(p) => p,
            Err(e) if e.kind() == io::ErrorKind::NotFound => {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("parent directory does not exist: {}", parent.display()),
                ).into());
            }
            Err(e) => return Err(e.into()),
        };
        let dir_name = self.output_dir.file_name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "output_dir ends with .. or has no final component",
            )
        })?;
        let output_dir = canon_parent.join(dir_name);

        let upper_cores = {
            let mut sys = System::new();
            sys.refresh_cpu_all();
            sys.cpus().len().saturating_sub(2).max(1)
        };

        let n_core = self.cores.unwrap_or(upper_cores).min(upper_cores);

        Ok(InitDedupBarcodeArgs::new(self.barcode_file, output_dir, self.tile_list, n_core))
    }
}

pub struct InitDedupBarcodeArgs {
    barcode_file: PathBuf,
    output_dir: PathBuf,
    tile_list: Vec<u64>,
    cores: usize,
}

impl InitDedupBarcodeArgs {
    #[inline]
    fn new(barcode_file: PathBuf, output_dir: PathBuf, tile_list: Vec<u64>, cores: usize) -> Self {
        Self {
            barcode_file,
            output_dir,
            tile_list,
            cores,
        }
    }

    #[inline]
    fn cores(&self) -> usize {
        self.cores
    }

    fn render_writer(&self) -> Result<(BufWriter<fs::File>, BufWriter<fs::File>), AppError> {
        // use for STAR to generate whitelist
        let barcode_whitelist = self.output_dir.join(format!("barcode_whitelist.txt"));
        let total_writer: BufWriter<fs::File> = BufWriter::new(
            fs::OpenOptions::new().create(true).write(true).open(barcode_whitelist)?
        );

        // use for map barcode to tile id
        let barcode_mapping = self.output_dir.join(format!("barcode_mapping.txt"));
        let map_writer = BufWriter::new(
            fs::OpenOptions::new().create(true).write(true).open(barcode_mapping)?
        );
        Ok((total_writer, map_writer))
    }

    fn parse_tile_record(line: &str) -> Result<(&str, &str, &str), AppError> {
        let mut f = line.splitn(4, '\t');
        f.next(); // chrom
        match (f.next(), f.next(), f.next()) {
            (Some(x), Some(y), Some(bc)) => Ok((bc, x, y)),
            _ =>
                Err(
                    AppError::Io(
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Invalid tile's barcode file format"
                        )
                    )
                ),
        }
    }

    pub fn dedup(self) -> Result<(), AppError> {
        fs::create_dir_all(&self.output_dir)?;

        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cores())
            .build()
            .expect("Failed to build Rayon thread pool");

        println!("Use {} threads for deduplicating", self.cores());

        // Set Writers
        let (mut total_writer, mut map_writer) = self.render_writer()?;

        // Set Step 2 ProgressBar
        let pb = ProgressBar::new(self.tile_list.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{prefix:.bold.dim} {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})"
                )
                .unwrap()
                .progress_chars("#>-")
        );
        pb.set_prefix("Processing tiles");

        let barcode_set = DashSet::new();
        let (sender, receiver) = crossbeam::channel::unbounded();

        let producer_handle = std::thread::spawn({
            let pb = pb.clone();
            move || {
                pool.install(|| {
                    self.tile_list.par_iter().try_for_each(|&tile_id| {
                        let mut reader = TbxReader::from_path(&self.barcode_file)?;
                        let tid = reader.tid(&tile_id.to_string())?;
                        reader.fetch(tid, 1000, 37100)?;

                        let mut writer = fs::OpenOptions
                            ::new()
                            .create(true)
                            .write(true)
                            .truncate(true)
                            .open(self.output_dir.join(format!("tile_{tile_id}.tsv.gz")))
                            .map(|f| GzEncoder::new(f, Compression::default()))
                            .map(BufWriter::new)?;
                        writeln!(writer, "cell_bc\tx_pos\ty_pos")?;

                        for record in reader.records() {
                            let record = unsafe { String::from_utf8_unchecked(record?) };

                            let (barcode, x, y) = Self::parse_tile_record(&record)?;

                            if barcode_set.insert(barcode.to_owned()) {
                                writeln!(writer, "{}\t{}\t{}", barcode, x, y)?;
                                sender
                                    .send((record.to_owned(), barcode.to_owned()))
                                    .map_err(|_| AppError::ChannelError)?;
                            }
                        }
                        pb.inc(1);
                        Ok::<(), AppError>(())
                    })
                })
            }
        });

        crossbeam
            ::scope(|s| {
                s.spawn(|_| {
                    for (record, barcode) in receiver {
                        writeln!(total_writer, "{}", barcode)?;
                        writeln!(map_writer, "{}", record)?;
                    }
                    Ok::<(), AppError>(())
                })
                    .join()
                    .unwrap()
            })
            .unwrap()?;

        producer_handle.join().unwrap()?;
        pb.finish_and_clear();

        Ok(())
    }
}
