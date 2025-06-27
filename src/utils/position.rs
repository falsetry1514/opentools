
use std::ops::Range;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Error, PartialEq)]
pub enum PositionError {
    #[error("Invalid format, expected 'read{{1/2}}:{{+/-}}:start-end'")]
    InvalidFormat,
    #[error("Invalid read specifier, must be 'read1' or 'read2'")]
    InvalidRead,
    #[error("Invalid strand, must be '+' or '-'")]
    InvalidStrand,
    #[error("Invalid start position, must be integer 0..150")]
    InvalidStart,
    #[error("Invalid end position, must be integer 0..150 or 'end'")]
    InvalidEnd,
    #[error("End position must be >= start position")]
    EndBeforeStart,
}

/// The struct stand for the position of sequence
#[derive(Debug, Copy, Clone)]
pub struct Position {
    /// false stand for read1, true stand for read2 
    read: bool,
    /// false stand for positive, true stand for negative
    strand: bool,
    /// Range in 0..150
    start: usize,
    /// Range in 0..150, must larger than start
    end: usize,
    /// The len of sequence
    len: usize
}

impl Position {
    pub fn new(read: bool, strand: bool, start: usize, end: usize) -> Self {
        let len = end - start;
        Self { read, strand, start, end, len }
    }

    #[inline]
    pub fn is_read2(&self) -> bool {self.read}

    #[inline]
    pub fn is_revcomp(&self) -> bool {self.strand}
    
    #[inline]
    pub fn start(&self) -> usize {self.start}

    #[inline]
    pub fn end(&self) -> usize {self.end}

    #[inline]
    pub fn len(&self) -> usize {self.len}

    #[inline]
    pub fn range(&self) -> Range<usize> {self.start..self.end}

    #[inline]
    pub fn safe_slice<'a, T>(&self, data: &'a [T]) -> &'a [T] {
        let start = std::cmp::min(self.start, data.len());
        let end = std::cmp::min(self.end, data.len());
        &data[start..end] // 自动处理越界
    }
}

impl FromStr for Position {
    type Err = PositionError;

    /// Parse the string from std::io::input into Position
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // split str into parts
        let parts: Vec<&str> = s.split(':').collect();
        if parts.len() != 3 {
            return Err(PositionError::InvalidFormat);
        }
        let range_parts: Vec<&str> = parts[2].split('-').collect();
        if range_parts.len() != 2 {
            return Err(PositionError::InvalidFormat);
        }

        // parse the part of read in position string ( read1 or read2 )
        let read = match parts[0] {
            "read1" => false,
            "read2" => true,
            _ => return Err(PositionError::InvalidRead),
        };

        // parse the part of strand in position string ( + or - )
        let strand = match parts[1] {
            "+" => false,
            "-" => true,
            _ => return Err(PositionError::InvalidStrand),
        };

        // parse the part of range in position string ( start-end )
        let start = match range_parts[0].parse::<usize>() {
            Err(_) => return Err(PositionError::InvalidStart),
            Ok(v) if v > 150 => return Err(PositionError::InvalidStart),
            Ok(v) => v,
        };
        let end = match range_parts[1].parse::<usize>() {
            Err(_) if range_parts[1].eq_ignore_ascii_case("end") => 150,
            Err(_) => return Err(PositionError::InvalidFormat),
            Ok(v) if v > 150 => return Err(PositionError::InvalidEnd),
            Ok(v) if v < start => return Err(PositionError::EndBeforeStart),
            Ok(v) => v,
        };

        Ok(Position::new(read, strand, start, end))
    }
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let read = if self.read { b'2' } else { b'1' };
        let strand = if self.strand { b'-' } else { b'+' };
        write!(f, "read{}:{}:{}-{}", read, strand, self.start, self.end)
    }
}
