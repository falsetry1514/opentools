
use crate::utils::{
    barcode_iter::{validate_absolute_filepath, validate_absolute_dirpath},
    error::AppError,
};
use crate::argparse::tilesmatch::is_valid_tile_id;
use std::fs;
use std::io::{self, Write, BufWriter};
use std::path::PathBuf;
use clap::Parser;
use dashmap::DashSet;
use rayon::prelude::*;
use rust_htslib::tbx::{self, Read};

#[derive(Parser, Debug)]
#[command(name = "dedupbarcode")]
pub struct DedupBarcodeArgs {
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
    tile_list: Vec<u64>,

    /// The path to the FASTQ file
    #[arg(
        short,
        long,
        required = true,
        value_parser = validate_absolute_dirpath,
    )]
    output_dir: PathBuf,
}

impl DedupBarcodeArgs {
    #[inline]
    pub fn tile_list(&self) -> &[u64] {
        &self.tile_list
    }

    pub fn dedup(self) -> Result<(), AppError> {
        let barcode_set = DashSet::new();

        // use for STAR to generate whitelist
        let barcode_whitelist = self.output_dir.join(format!("barcode_whitelist.txt"));
        let mut total_writer = BufWriter::new(
            fs::OpenOptions::new().create(true).write(true).open(barcode_whitelist)?
        );

        // use for map barcode to tile id
        let barcode_mapping = self.output_dir.join(format!("barcode_mapping.txt"));
        let mut map_writer = BufWriter::new(
            fs::OpenOptions::new().create(true).write(true).open(barcode_mapping)?
        );

        let (sender, receiver) = crossbeam::channel::unbounded();
    
        let producer_handle = std::thread::spawn(
            move || {
                self.tile_list.par_iter().try_for_each(|&tile_id| {
                    let tile_file = self.output_dir.join(format!("{tile_id}.txt"));
                    let mut writer = BufWriter::new(
                        fs::OpenOptions::new().create(true).write(true).open(tile_file)?
                    );
        
                    let mut reader = tbx::Reader::from_path(&self.barcode_file)?;
                    let tid = reader.tid(&tile_id.to_string())?;
                    reader.fetch(tid, 1000, 37100)?;

                    writeln!(writer, "tile_id\tx_po\ty_pos\tbarcode")?;
                    for record in reader.records() {
                        let record = record?;
                        let record = unsafe { String::from_utf8_unchecked(record) };
                        let barcode = record.splitn(4, '\t').nth(3).ok_or(AppError::IoError(
                            io::Error::new(io::ErrorKind::InvalidData, "Invalid tile's barcode file format")
                        ))?;

                        if barcode_set.insert(barcode.to_string()) {
                            writeln!(writer, "{}", record)?;
                            sender.send((record.to_owned(), barcode.to_string())).map_err(|_| AppError::ChannelError)?;
                        }
                    }
                    Ok::<(), AppError>(())
                })
            }
        );

        crossbeam::scope(|s| {
            s.spawn(|_| {
                for (record, barcode) in receiver {
                    writeln!(total_writer, "{}", barcode)?;
                    writeln!(map_writer, "{}", record)?;
                }
                Ok::<(), AppError>(())
            }).join().unwrap()
        }).unwrap()?;

        producer_handle.join().unwrap()?;
        
        Ok(())
    }
}