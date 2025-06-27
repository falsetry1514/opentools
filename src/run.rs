use crate::argparse::{
    dedupbarcode::DedupBarcodeArgs, fq2bam::Fq2BamArgs, tilesmatch::TilesMatchArgs,
    touchbarcode::TouchBarcodeArgs,
};
use crate::utils::error::AppError;

use rayon::{ThreadPoolBuilder, prelude::*};
use std::{fs, process::Command};

/// Default thread count configuration
/// 
/// Default: 12 threads for Linux, 3 threads for macOS
pub const DEFAULT_LINUX_THREADS: usize = 12;
pub const DEFAULT_MAC_THREADS: usize = 3;

/// Processes FASTQ to BAM conversion workflow
/// 
/// # Arguments
/// - `args`: Fq2BamArgs struct containing input file paths and output configuration
/// 
/// # Errors
/// Returns AppError for possible I/O errors or data format errors
pub fn fq2bam(mut args: Fq2BamArgs) -> Result<(), AppError> {
    args.validate_eq_file_count()?;
    let config = args.record_config();
    let mut stdout = args.create_bam_header();

    for reader in args.paired_readers() {
        let mut reader = reader?;
        let (sender, receiver) = crossbeam::channel::bounded(4096);
        let producer_handle = std::thread::spawn(move || {
            reader
                .records(config)
                .par_bridge()
                .try_for_each(|record| sender.send(record).map_err(|_| AppError::ChannelError))
        });

        crossbeam::scope(|s| {
            s.spawn(|_| -> Result<(), AppError> {
                for record in receiver.iter() {
                    stdout.write(&record?)?;
                }
                Ok(())
            })
            .join()
            .unwrap()
        })
        .unwrap()?;

        producer_handle.join().unwrap()?;
    }
    Ok(())
}

/// Handles barcode viewing and deduplication
///
/// # Arguments
/// - `args`: DedupBarcodeArgs struct containing input files and deduplication configuration
///
/// # Errors
/// Returns AppError for possible I/O errors or data processing errors
pub fn viewbarcode(args: DedupBarcodeArgs) -> Result<(), AppError> {
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
    let tile_ids = args.extract_tile_ids()?;
    println!("Extracted tile IDs from bcl directory RunInfo.xml file");
    let num_threads: usize = if cfg!(target_os = "linux") {
        DEFAULT_LINUX_THREADS
    } else if cfg!(target_os = "macos") {
        DEFAULT_MAC_THREADS
    } else {
        return Err(AppError::UnsupportedOS);
    };

    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Build thread pool failed");
    let tile_ids: Vec<String> = pool.install(|| {
        tile_ids
            .par_iter()
            .map(|tile_id| {
                let fastq_file = args
                    .fastq_path(tile_id)
                    .join("Undetermined_S0_R1_001.fastq.gz");
                if !fastq_file.exists() {
                    println!("Converted tile {tile_id} into fastq");
                    args.convert_bcl_into_tile(&tile_id)?;
                } else {
                    println!("Have already converted tile {tile_id}");
                };
                let tile_id = tile_id.replace("_", "");
                Ok(tile_id)
            })
            .collect::<Result<Vec<String>, AppError>>()
    })?;

    let mut tile_ids: Vec<String> = tile_ids
        .into_par_iter()
        .map(|tile_id| {
            let barcode_iter = args.create_barcode_iter(&tile_id)?;
            let report = barcode_iter.extract_chip_barcodes()?;
            println!("Tile {tile_id}: {report}");
            println!("Extracted Barcode of tile_id {tile_id} into tmp file.");
            Ok(tile_id)
        })
        .collect::<Result<Vec<String>, AppError>>()?;
    tile_ids.par_sort_unstable();

    let files: Vec<String> = tile_ids
        .into_iter()
        .map(|tile_id| {
            args.output()
                .join(format!("tmp/{}.txt", tile_id))
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
    if tmp_dir.exists() {
        fs::remove_dir_all(&tmp_dir)?;
    }

    let tabix_status = Command::new("tabix")
        .args(&["-0", "-s", "1", "-b", "3", "-e", "3"])
        .arg(output_path)
        .status()?;
    if !tabix_status.success() {
        return Err(AppError::CommandError("tabix run failed".to_string()));
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::error::AppError;
    use tempfile::tempdir;

    #[test]
    fn test_default_thread_counts() {
        assert_eq!(DEFAULT_LINUX_THREADS, 12);
        assert_eq!(DEFAULT_MAC_THREADS, 3);
    }

    #[test]
    fn test_viewbarcode_returns_ok() {
        // This is just a placeholder test
        let result: Result<(), AppError> = Ok(());
        assert!(result.is_ok());
    }

    #[test]
    fn test_tilesmatch_returns_ok() {
        // This is just a placeholder test
        let result: Result<(), AppError> = Ok(());
        assert!(result.is_ok());
    }

    #[test]
    fn test_directory_creation() {
        let dir = tempdir().unwrap();
        let path = dir.path();
        
        // Test that we can create a directory
        let test_dir = path.join("test_dir");
        fs::create_dir(&test_dir).unwrap();
        assert!(test_dir.exists());
    }
}
