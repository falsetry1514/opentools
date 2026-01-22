use clap::Parser;
use opentools::argparse::{ Cli, Commands };
use opentools::run;
use opentools::utils::error::AppError;

fn main() -> Result<(), AppError> {
    let cli = Cli::parse();

    match cli.command {
        Commands::TouchBarcode(args) => run::touchbarcode(args)?,
        Commands::ViewBarcode(args) => run::dedupbarcode(args)?,
        Commands::TilesMatch(args) => run::tilesmatch(args)?,
        Commands::DemuxI2(args) => run::demux_i2(args)?,
    }

    Ok(())
}
