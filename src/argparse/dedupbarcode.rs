// ------------------
// 1. std
// ------------------
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

// ------------------
// 2. external crates
// ------------------
use clap::{Parser, ValueEnum};
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use rust_htslib::tbx::{Read, Reader as TbxReader};
use sysinfo::System;

// ------------------
// 3. inner crate
// ------------------
use crate::utils::{
    barcode_iter::validate_absolute_filepath, error::AppError, tile_ids::is_valid_tile_id,
};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum AmbiguityPolicy {
    /// Keep every record, even when one barcode maps to multiple distinct spatial coordinates.
    KeepAll,
    /// Drop every record of a barcode that maps to multiple distinct spatial coordinates.
    DropAll,
}

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

    /// How to handle a barcode that maps to multiple distinct spatial coordinates
    #[arg(long, value_enum, default_value = "keep-all")]
    ambiguity_policy: AmbiguityPolicy,

    /// Emit per-tile tsv.gz files
    #[arg(long)]
    tile_files: bool,

    /// Do not emit barcode_whitelist.txt / barcode_mapping.txt summary files
    #[arg(long)]
    no_summary_files: bool,
}

impl DedupBarcodeArgs {
    pub fn init(self) -> Result<InitDedupBarcodeArgs, AppError> {
        if !self.tile_files && self.no_summary_files {
            return Err(AppError::Io(io::Error::new(
                io::ErrorKind::InvalidInput,
                "at least one output must be enabled: pass --tile-files and/or drop --no-summary-files",
            )));
        }

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
                )
                .into());
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

        Ok(InitDedupBarcodeArgs::new(
            self.barcode_file,
            output_dir,
            self.tile_list,
            n_core,
            self.ambiguity_policy,
            self.tile_files,
            !self.no_summary_files,
        ))
    }
}

pub struct InitDedupBarcodeArgs {
    barcode_file: PathBuf,
    output_dir: PathBuf,
    tile_list: Vec<u64>,
    cores: usize,
    ambiguity_policy: AmbiguityPolicy,
    emit_tile_files: bool,
    emit_summary_files: bool,
}

impl InitDedupBarcodeArgs {
    #[inline]
    fn new(
        barcode_file: PathBuf,
        output_dir: PathBuf,
        tile_list: Vec<u64>,
        cores: usize,
        ambiguity_policy: AmbiguityPolicy,
        emit_tile_files: bool,
        emit_summary_files: bool,
    ) -> Self {
        Self {
            barcode_file,
            output_dir,
            tile_list,
            cores,
            ambiguity_policy,
            emit_tile_files,
            emit_summary_files,
        }
    }

    #[inline]
    fn cores(&self) -> usize {
        self.cores
    }

    fn parse_tile_record(line: &str) -> Result<(&str, &str, &str), AppError> {
        let mut f = line.splitn(4, '\t');
        f.next(); // chrom
        match (f.next(), f.next(), f.next()) {
            (Some(x), Some(y), Some(bc)) => Ok((bc, x, y)),
            _ => Err(AppError::Io(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid tile's barcode file format",
            ))),
        }
    }

    fn scan_ambiguous(&self, pool: &rayon::ThreadPool) -> Result<HashSet<String>, AppError> {
        let barcode_file = &self.barcode_file;
        let tile_list = &self.tile_list;

        type Sight = Option<(u64, String, String)>;
        type LocalMap = HashMap<String, Sight>;

        // barcode -> Some(first coordinate) | None (= a different coordinate was observed)
        let map = pool.install(|| -> Result<LocalMap, AppError> {
            tile_list
                .par_iter()
                .try_fold(
                    LocalMap::new,
                    |mut local, &tile_id| -> Result<_, AppError> {
                        let mut reader = TbxReader::from_path(barcode_file)?;
                        let tid = reader.tid(&tile_id.to_string())?;
                        reader.fetch(tid, 1000, 37100)?;
                        for record in reader.records() {
                            let record = unsafe { String::from_utf8_unchecked(record?) };
                            let (bc, x, y) = Self::parse_tile_record(&record)?;
                            match local.get_mut(bc) {
                                Some(Some((ft, fx, fy)))
                                    if *ft == tile_id && fx == x && fy == y => {}
                                Some(slot) => {
                                    *slot = None;
                                }
                                None => {
                                    local.insert(
                                        bc.to_owned(),
                                        Some((tile_id, x.to_owned(), y.to_owned())),
                                    );
                                }
                            }
                        }
                        Ok(local)
                    },
                )
                .try_reduce(LocalMap::new, |mut a, b| {
                    for (k, v) in b {
                        match a.entry(k) {
                            Entry::Occupied(mut e) => {
                                e.insert(None);
                            }
                            Entry::Vacant(e) => {
                                e.insert(v);
                            }
                        }
                    }
                    Ok(a)
                })
        })?;

        let amb: HashSet<String> = map
            .into_iter()
            .filter_map(|(k, v)| if v.is_none() { Some(k) } else { None })
            .collect();
        Ok(amb)
    }

    pub fn dedup(self) -> Result<(), AppError> {
        if !self.emit_tile_files && !self.emit_summary_files {
            return Err(
                AppError::Io(
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "at least one output must be enabled: --tile-files and/or --no-summary-files=false"
                    )
                )
            );
        }

        fs::create_dir_all(&self.output_dir)?;

        let pool = ThreadPoolBuilder::new()
            .num_threads(self.cores())
            .build()
            .expect("Failed to build Rayon thread pool");

        println!("Use {} threads for deduplicating", self.cores());
        println!("Ambiguity policy: {:?}", self.ambiguity_policy);
        println!("Emit per-tile files: {}", self.emit_tile_files);
        println!("Emit summary files: {}", self.emit_summary_files);

        // Pass 1: only built when DropAll is requested
        let ambiguous: Option<Arc<HashSet<String>>> = match self.ambiguity_policy {
            AmbiguityPolicy::KeepAll => None,
            AmbiguityPolicy::DropAll => {
                println!("Scanning for barcodes with multiple spatial coordinates...");
                let amb = self.scan_ambiguous(&pool)?;
                println!("Found {} ambiguous barcode(s) to drop", amb.len());
                Some(Arc::new(amb))
            }
        };

        let (mut total_writer, mut map_writer) = if self.emit_summary_files {
            let bw = BufWriter::new(
                fs::OpenOptions
                    ::new()
                    .create(true)
                    .write(true)
                    .truncate(true)
                    .open(self.output_dir.join("barcode_whitelist.txt"))?
            );
            let mw = BufWriter::new(
                fs::OpenOptions
                    ::new()
                    .create(true)
                    .write(true)
                    .truncate(true)
                    .open(self.output_dir.join("barcode_mapping.txt"))?
            );
            (Some(bw), Some(mw))
        } else {
            (None, None)
        };

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

        let (sender, receiver) = crossbeam::channel::unbounded();

        let emit_tile_files = self.emit_tile_files;
        let emit_summary_files = self.emit_summary_files;
        let output_dir = self.output_dir.clone();

        let producer_handle = std::thread::spawn({
            let pb = pb.clone();
            let ambiguous = ambiguous.clone();
            move || {
                pool.install(|| {
                    self.tile_list.par_iter().try_for_each(|&tile_id| {
                        let mut reader = TbxReader::from_path(&self.barcode_file)?;
                        let tid = reader.tid(&tile_id.to_string())?;
                        reader.fetch(tid, 1000, 37100)?;

                        let mut writer = if emit_tile_files {
                            let f = fs::OpenOptions
                                ::new()
                                .create(true)
                                .write(true)
                                .truncate(true)
                                .open(output_dir.join(format!("tile_{tile_id}.tsv.gz")))?;
                            let mut w = BufWriter::new(GzEncoder::new(f, Compression::default()));
                            writeln!(w, "cell_bc\tx_pos\ty_pos")?;
                            Some(w)
                        } else {
                            None
                        };

                        for record in reader.records() {
                            let record = unsafe { String::from_utf8_unchecked(record?) };
                            let (barcode, x, y) = Self::parse_tile_record(&record)?;

                            // Under DropAll, skip every record of an ambiguous barcode
                            if let Some(amb) = ambiguous.as_ref() {
                                if amb.contains(barcode) {
                                    continue;
                                }
                            }

                            if let Some(w) = writer.as_mut() {
                                writeln!(w, "{}\t{}\t{}", barcode, x, y)?;
                            }
                            if emit_summary_files {
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

        if emit_summary_files {
            let mut total_writer = total_writer.take().unwrap();
            let mut map_writer = map_writer.take().unwrap();
            crossbeam
                ::scope(|s| {
                    s.spawn(move |_| {
                        let mut seen_total: HashSet<String> = HashSet::new();
                        for (record, barcode) in receiver {
                            // Whitelist is deduplicated by barcode (globally)
                            if seen_total.insert(barcode.clone()) {
                                writeln!(total_writer, "{}", barcode)?;
                            }
                            // Mapping keeps every record that survived the ambiguity filter
                            writeln!(map_writer, "{}", record)?;
                        }
                        Ok::<(), AppError>(())
                    })
                        .join()
                        .unwrap()
                })
                .unwrap()?;
        } else {
            // No summary output: producer never sends, so just drop the receiver
            drop(receiver);
        }

        producer_handle.join().unwrap()?;
        pb.finish_and_clear();

        Ok(())
    }
}
