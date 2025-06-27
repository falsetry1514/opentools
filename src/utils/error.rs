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
    
    /// Mismatched number of input files: {n1} forward files, {n2} reverse files
    #[error("Mismatched number of input files: {n1} forward files, {n2} reverse files")]
    NotEqualFileNumber { n1: usize, n2: usize },
    
    /// Inconsistent line counts in Fastq files: file1 at line {line1}, file2 at line {line2}
    #[error("Inconsistent line counts in Fastq files: file1 at line {line1}, file2 at line {line2}")]
    FastqFileLengthNotEqual { line1: u64, line2: u64 },
    
    /// Fastq parsing error: {0}
    #[error("Fastq parsing error: {0}")]
    FastqParseError(#[source] SeqIoError),
    
    /// Mismatched Fastq IDs at line {line}: ID1={id1}, ID2={id2}
    #[error("Mismatched Fastq IDs at line {line}: ID1={id1}, ID2={id2}")]
    FastqIdMismatch { line: u64, id1: String, id2: String },
    
    /// BAM record operation error: {0}
    #[error("BAM record operation error: {0}")]
    BamRecordError(#[from] BamError),
    
    /// Empty tile IDs list: {0:?}
    #[error("Empty tile IDs list: {0:?}")]
    EmptyTileIDsList(PathBuf),
    
    /// Invalid barcode pattern: {0}
    #[error("Invalid barcode pattern: {0}")]
    InvalidBarcodePattern(String),
    
    /// Barcode contains invalid UTF-8 characters
    #[error("Barcode contains invalid UTF-8 characters")]
    InvalidUtf8InBarcode,
    
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    #[test]
    fn test_io_error_display() {
        let io_error = io::Error::new(io::ErrorKind::NotFound, "file not found");
        let app_error = AppError::IoError(io_error);
        assert_eq!(
            app_error.to_string(),
            "IO operation error: file not found"
        );
    }

    #[test]
    fn test_not_equal_file_number_display() {
        let app_error = AppError::NotEqualFileNumber { n1: 3, n2: 5 };
        assert_eq!(
            app_error.to_string(),
            "Mismatched number of input files: 3 forward files, 5 reverse files"
        );
    }

    #[test]
    fn test_command_not_found_display() {
        let app_error = AppError::CommandNotFound("my_command".to_string());
        assert_eq!(
            app_error.to_string(),
            "System command not found: my_command"
        );
    }
}
