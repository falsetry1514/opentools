use std::path::PathBuf;
use thiserror::Error;
use seq_io::fastq::{ Error as SeqIoError, ErrorPosition };

#[derive(Debug, Error)]
pub enum FastqError {
    /// IO error
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    /// Sequence and quality lengths are not equal
    #[error("Unequal sequence/quality lengths: seq={seq}, qual={qual}, pos={pos:?}")]
    UnequalLengths {
        seq: usize,
        qual: usize,
        pos: ErrorPosition,
    },

    /// Invalid start byte (expected '@')
    #[error("Invalid start byte {found} at pos {pos:?}")]
    InvalidStart {
        found: u8,
        pos: ErrorPosition,
    },

    /// Invalid separator byte (expected '+')
    #[error("Invalid separator byte {found} at pos {pos:?}")]
    InvalidSep {
        found: u8,
        pos: ErrorPosition,
    },

    /// Truncated record found
    #[error("Unexpected end of file at pos {pos:?}")]
    UnexpectedEnd {
        pos: ErrorPosition,
    },

    /// Buffer size limit reached
    #[error("Buffer size limit reached")]
    BufferLimit,
}

/// Unified error handling type for the application
///
/// Uses thiserror for deriving error handling, providing clear error context information
#[derive(Debug, Error)]
pub enum AppError {
    /// IO operation error: {0}
    #[error("IO operation error: {0}")]
    Io(#[from] std::io::Error),

    /// FASTQ parsing error
    #[error("FASTQ error: {0}")]
    Fastq(#[from] FastqError),

    /// BAM record operation error: {0}
    #[error("BAM error: {0}")]
    Bam(#[from] rust_htslib::errors::Error),

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

    /// Early stop
    #[error("Early stop")]
    EarlyStop,
}

impl From<SeqIoError> for FastqError {
    fn from(err: SeqIoError) -> Self {
        match err {
            SeqIoError::Io(err) => FastqError::Io(err),
            SeqIoError::UnequalLengths { seq, qual, pos } =>
                FastqError::UnequalLengths { seq, qual, pos },
            SeqIoError::InvalidStart { found, pos } => FastqError::InvalidStart { found, pos },
            SeqIoError::InvalidSep { found, pos } => FastqError::InvalidSep { found, pos },
            SeqIoError::UnexpectedEnd { pos } => FastqError::UnexpectedEnd { pos },
            SeqIoError::BufferLimit => FastqError::BufferLimit,
            // _ => AppError::FastqParseError(err),
        }
    }
}

impl From<SeqIoError> for AppError {
    fn from(err: SeqIoError) -> Self {
        AppError::Fastq(FastqError::from(err))
    }
}
