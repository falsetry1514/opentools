use std::path::PathBuf;
// ------------------
// 1. std
// ------------------
use std::sync::{ Arc, Mutex };
use std::fs;
use std::process::Command;
use std::io::{ Write, BufReader };

// ------------------
// 2. external crates
// ------------------
use rust_htslib::bgzf;
use rust_htslib::tpool;
use indicatif::{ ProgressBar, ProgressStyle };
use rayon::{ ThreadPoolBuilder, prelude::* };

// ------------------
// 3. inner crate
// ------------------
use crate::argparse::{
    dedupbarcode::DedupBarcodeArgs,
    tilesmatch::TilesMatchArgs,
    touchbarcode::TouchBarcodeArgs,
    add_bc::AddBCArgs,
    stats_bc::StatsBCArgs,
};
use crate::utils::error::AppError;

/// Handles barcode viewing and deduplication
///
/// # Arguments
/// - `args`: DedupBarcodeArgs struct containing input files and deduplication configuration
///
/// # Errors
/// Returns AppError for possible I/O errors or data processing errors
pub fn dedupbarcode(args: DedupBarcodeArgs) -> Result<(), AppError> {
    let args = args.init()?;
    args.dedup()?;
    Ok(())
}

/// Handles barcode preprocessing workflow
///
/// # Arguments
/// - `args`: TouchBarcodeArgs struct containing input path and output configuration
///
/// # Errors
/// Returns AppError for possible I/O errors, system command not found, or execution failure
pub fn touchbarcode(args: TouchBarcodeArgs) -> Result<(), AppError> {
    let args = args.init()?;
    args.validate_command()?;

    // Create output directories
    let fastq_dir = args.output().join("fastq");
    let tmp_dir = args.output().join("tmp");
    if !fastq_dir.exists() {
        fs::create_dir(&fastq_dir)?;
    }
    if !tmp_dir.exists() {
        fs::create_dir(&tmp_dir)?;
    }

    // Extract tile IDs
    let mut tile_ids = args.extract_tile_ids()?;
    println!("Extracted tile IDs from bcl directory RunInfo.xml file");

    // Set cores
    println!("Use {} threads for conversion", args.cores());

    let pool = ThreadPoolBuilder::new()
        .num_threads(args.cores())
        .build()
        .expect("Build thread pool failed");

    // Set Step 1 ProgressBar
    let pb = ProgressBar::new(tile_ids.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{prefix:.bold.dim} {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})"
            )
            .unwrap()
            .progress_chars("#>-")
    );
    pb.set_prefix("BCL→FASTQ Conversion");

    // Step 1
    pool.install(|| {
        tile_ids
            .par_iter()
            .map(|tile_id| {
                let fastq_file = args.fastq_path(tile_id).join("Undetermined_S0_R1_001.fastq.gz");

                if !fastq_file.exists() {
                    args.convert_bcl_into_tile(tile_id)?;
                }

                pb.inc(1);
                Ok(())
            })
            .collect::<Result<(), AppError>>()
    })?;
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.bold.dim} [{elapsed_precise}] {msg}")
            .unwrap()
    );
    pb.finish_with_message("Step 1: All tiles processed! ✅");

    // Set Step 2 ProgressBar
    let pb = ProgressBar::new(tile_ids.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{prefix:.bold.dim} {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} {msg} ({eta})"
            )
            .unwrap()
            .progress_chars("#>-")
    );
    pb.set_prefix("Extracting barcodes");

    // Set Step 2 logfile
    let log_path = args.output().join("touchbarcode.log");
    let log_file = fs::OpenOptions
        ::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&log_path)
        .map(Mutex::new)
        .map(Arc::new)?;

    // Step 2
    pool.install(|| {
        tile_ids
            .par_iter()
            .map(|tile_id| {
                let barcode_iter = args.create_barcode_iter(tile_id)?;
                let report = barcode_iter.extract_chip_barcodes(args.render_writer(tile_id)?)?;

                {
                    let mut file = log_file.lock().unwrap();
                    writeln!(file, "Tile {tile_id}: {report}")?;
                }

                pb.inc(1);
                Ok(())
            })
            .collect::<Result<(), AppError>>()
    })?;
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.bold.dim} [{elapsed_precise}] {msg}")
            .unwrap()
    );
    pb.finish_with_message("Step 2: Barcode extraction completed! ✅");

    // Step 2.5: sort
    tile_ids.par_sort_unstable();
    let files: Vec<PathBuf> = tile_ids
        .iter()
        .map(|tile_id| { args.output().join(format!("tmp/{tile_id}.txt")) })
        .collect();

    // Step 3: merge
    let output_path = args.output().join("barcodes.txt.gz");
    let mut writer = bgzf::Writer::from_path_with_level(&output_path, bgzf::CompressionLevel::Maximum)?;
    let htslib_pool = tpool::ThreadPool::new(args.cores() as u32)?;
    writer.set_thread_pool(&htslib_pool)?;
    writer.write_all(b"#tile_id\tx_pos\ty_pos\tbarcode\n")?;

    let pb = ProgressBar::new(files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{prefix:.bold.dim} {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} {msg} ({eta})"
            )
            .unwrap()
            .progress_chars("#>-")
    );
    pb.set_prefix("Merging tiles");

    for path in &files {
        let mut reader = fs::File::open(path).map(BufReader::new)?;
        std::io::copy(&mut reader, &mut writer)?;
        pb.inc(1);
    }
    writer.flush()?;
    drop(writer);

    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.bold.dim} [{elapsed_precise}] {msg}")
            .unwrap()
    );
    pb.finish_with_message("Step 3: Merging tiles completed! ✅");

    println!("Barcodes written to: {}", output_path.display());
    if tmp_dir.exists() {
        fs::remove_dir_all(&tmp_dir)?;
    }

    // Step 4: write tabix index
    let tabix_status = Command::new("tabix")
        .args(&["-0", "-s", "1", "-b", "3", "-e", "3", "-@", args.cores().to_string().as_ref()])
        .arg(&output_path)
        .status()?;
    if !tabix_status.success() {
        return Err(AppError::CommandError("tabix run failed".to_string()));
    }
    println!("Created index for: {}", output_path.display());
    Ok(())
}

/// Handles tile matching analysis
///
/// # Arguments
/// - `args`: TilesMatchArgs struct containing matching threshold and input files
///
/// # Errors
/// Returns AppError for possible I/O errors or data processing errors
pub fn tilesmatch(args: TilesMatchArgs) -> Result<(), AppError> {
    let args = args.init()?;

    let reports = args.search_tile()?;

    if !args.quiet() {
        println!("Tile id\tTotal number\tMatched number\tMatch ratio\tPass threshold");
    }

    reports.into_iter().for_each(|report| {
        if args.quiet() {
            if report.pass_threshold() {
                print!("{} ", report.tile_id());
            }
        } else {
            println!("{report}")
        }
    });

    Ok(())
}

pub fn add_bc(args: AddBCArgs) -> Result<(), AppError> {
    let args = args.init()?;
    args.run_add_cb()?;
    Ok(())
}

pub fn stats_bc(args: StatsBCArgs) -> Result<(), AppError> {
    let args = args.init()?;
    args.run_stats_cb()?;
    Ok(())
}
