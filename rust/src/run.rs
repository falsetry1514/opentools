use crate::argparse::{
    dedupbarcode::DedupBarcodeArgs, 
    tilesmatch::TilesMatchArgs,
    touchbarcode::TouchBarcodeArgs,
};
use crate::utils::error::AppError;

use indicatif::{ProgressBar, ProgressStyle};
use rayon::{ThreadPoolBuilder, prelude::*};
use sysinfo::System;
use std::sync::{Arc, Mutex};
use std::{fs, process::Command};
use std::io::Write;

/// Handles barcode viewing and deduplication
///
/// # Arguments
/// - `args`: DedupBarcodeArgs struct containing input files and deduplication configuration
///
/// # Errors
/// Returns AppError for possible I/O errors or data processing errors
pub fn dedupbarcode(args: DedupBarcodeArgs) -> Result<(), AppError> {
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
    let args = args.init();
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
    let num_threads: usize = {
        let mut sys = System::new();

        sys.refresh_memory();
        sys.refresh_cpu_usage();

        let avail_gb = sys.available_memory() / 1024 / 1024 / 1024;
        let usable = (avail_gb * 4 / 5 / 5) as usize;
        let cores = sys.cpus().len().max(1);

        usable.min(cores)
    };
    println!("Use {} threads for conversion", num_threads);

    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Build thread pool failed");

    let pb = ProgressBar::new(tile_ids.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.bold.dim} {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    pb.set_prefix("BCL→FASTQ Conversion");

    pool.install(|| {
        tile_ids
            .par_iter()
            .map(|tile_id| {
                let fastq_file = args
                    .fastq_path(tile_id)
                    .join("Undetermined_S0_R1_001.fastq.gz");

                if !fastq_file.exists() {
                    args.convert_bcl_into_tile(tile_id)?;
                }

                pb.inc(1);
                Ok(())
            })
            .collect::<Result<(), AppError>>()
    })?;
    pb.finish_with_message("Step 1: All tiles processed! ✅");

    let pb = ProgressBar::new(tile_ids.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{prefix:.bold.dim} {spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );
    pb.set_prefix("Extracting barcodes");

    let log_path = args.output().join("touchbarcode.log");
    let log_file = Arc::new(Mutex::new(
        fs::OpenOptions::new().create(true).append(true).open(&log_path)?
    ));

    tile_ids
        .par_iter()
        .map(|tile_id| {
            let barcode_iter = args.create_barcode_iter(tile_id)?;
            let report = barcode_iter.extract_chip_barcodes()?;

            {
                let mut file = log_file.lock().unwrap();
                writeln!(file, "Tile {tile_id}: {report}")?;
            }

            pb.inc(1);
            Ok(())
        })
        .collect::<Result<(), AppError>>()?;
    pb.finish_with_message("Step 2: Barcode extraction completed! ✅");
    tile_ids.par_sort_unstable();

    let files: Vec<String> = tile_ids
        .iter()
        .map(|tile_id| {
            args.output()
                .join(format!("tmp/{tile_id}.txt"))
                .display()
                .to_string()
        })
        .collect();
    let output_path = args.output().join("barcodes.txt.gz");

    let output = Command::new("bash")
        .arg("-c")
        .arg(&format!(
            "{{ echo '#tile_id\tx_pos\ty_pos\tbarcode'; cat {}; }} | bgzip -@ $(nproc) > {}",
            files.join(" "),
            output_path.display()
        ))
        .output()?;
    if !output.status.success() {
        return Err(AppError::CommandError(format!(
            "bgzip run failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    println!("Barcodes written to: {}", output_path.display());
    if tmp_dir.exists() {
        fs::remove_dir_all(&tmp_dir)?;
    }

    let tabix_status = Command::new("tabix")
        .args(&["-0", "-s", "1", "-b", "3", "-e", "3", "-@", "$(nproc)"])
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
        println!("Tile id\tTotal number\tMatched number\tMatch ratio\tPass threshold")
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
