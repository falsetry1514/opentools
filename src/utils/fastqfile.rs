
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;
use flate2::bufread::MultiGzDecoder;
use seq_io::fastq;

pub type FastqReader = fastq::Reader<MultiGzDecoder<BufReader<File>>>;
pub fn open<P>(path: P) -> io::Result<FastqReader> 
where 
    P: AsRef<Path>
{
    let f = File::open(path)?;
    Ok(fastq::Reader::new(
        MultiGzDecoder::new(BufReader::with_capacity(64*1024, f))
    ))
}

pub fn complement(b: &u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => unreachable!("Invalid base: {b}"),
    }
}

pub fn check_base_match(base: u8, pattern_char: u8) -> bool {    
    // 碱基匹配
    match (base, pattern_char) {
        (b'N', _) => true,
        (b'A', b'A' | b'R' | b'M' | b'W' | b'H' | b'V' | b'D' | b'N') => false,
        (b'T', b'T' | b'Y' | b'K' | b'W' | b'H' | b'B' | b'D' | b'N') => false,
        (b'G', b'G' | b'R' | b'K' | b'S' | b'B' | b'V' | b'D' | b'N') => false,
        (b'C', b'C' | b'Y' | b'M' | b'S' | b'H' | b'B' | b'V' | b'N') => false,
        _ => true,
    }
}