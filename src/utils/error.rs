use std::path::PathBuf;
use thiserror::Error;
use seq_io::fastq::Error as SeqIoError;
use rust_htslib::errors::Error as BamError;

/// Unified error handling type for the application
/// 
/// Uses thiserror for deriving error handling, providing clear error context information
#[derive(Debug, Error)]
pub enum AppError {
    /// IO operation error: {0}
    #[error("IO operation error: {0}")]
    IoError(#[from] std::io::Error),
    
    /// Fastq parsing error: {0}
    #[error("Fastq parsing error: {0}")]
    FastqParseError(#[source] SeqIoError),
    
    /// BAM record operation error: {0}
    #[error("BAM record operation error: {0}")]
    BamRecordError(#[from] BamError),
    
    /// Empty tile IDs list: {0:?}
    #[error("Empty tile IDs list: {0:?}")]
    EmptyTileIDsList(PathBuf),
    
    /// Invalid barcode pattern: {0}
    #[error("Invalid barcode pattern: {0}")]
    InvalidBarcodePattern(String),
    
    /// Thread channel communication failed
    #[error("Thread channel communication failed")]
    ChannelError,
    
    /// Unsupported operating system
    #[error("Unsupported operating system")]
    UnsupportedOS,
    
    /// Docker image not found: {0}
    #[error("Docker image not found: {0}")]
    DockerImageNotFound(String),
    
    /// System command not found: {0}
    #[error("System command not found: {0}")]
    CommandNotFound(String),
    
    /// Command execution failed: {0}
    #[error("Command execution failed: {0}")]
    CommandError(String),
}

impl From<SeqIoError> for AppError {
    fn from(err: SeqIoError) -> Self {
        match err {
            SeqIoError::Io(err) => AppError::IoError(err),
            _ => AppError::FastqParseError(err),
        }
    }
}