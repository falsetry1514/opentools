
pub mod touchbarcode;
pub mod dedupbarcode;
pub mod tilesmatch;

use clap::{Parser, Subcommand};
use self::{
    touchbarcode::TouchBarcodeArgs,
    dedupbarcode::DedupBarcodeArgs,
    tilesmatch::TilesMatchArgs,
};

/// Command line arguments resolve the main structure
/// 
/// Use the clap-derived macro to implement command line parameter parsing
#[derive(Parser)]
#[command(name = "opentools")]
#[command(version = "1.0")]
#[command(about = "OpenST toolbox", long_about = None)]
#[command(next_line_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Subcommand enumeration definitions
/// 
/// Each variant corresponds to a specific tool function
#[derive(Subcommand)]
pub enum Commands {
    #[clap(name="touchbarcode")]
    TouchBarcode(TouchBarcodeArgs),
    #[clap(name="dedupbarcode")]
    ViewBarcode(DedupBarcodeArgs),
    #[clap(name="tilesmatch")]
    TilesMatch(TilesMatchArgs),
}
