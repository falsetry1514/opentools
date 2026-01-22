use std::fs::File;
use std::io::{ self, BufReader, Read };
use std::path::Path;
use flate2::bufread::MultiGzDecoder;
use seq_io::fastq::Reader as FastqReader;

// pub type PlainFastqReader = FastqReader<BufReader<File>>;
// pub fn plain_open<P>(path: P) -> io::Result<PlainFastqReader> where P: AsRef<Path> {
//     let f = File::open(path)?;
//     Ok(FastqReader::new(BufReader::with_capacity(64 * 1024, f)))
// }

// pub type GzFastqReader = FastqReader<MultiGzDecoder<BufReader<File>>>;
// pub fn gz_open<P>(path: P) -> io::Result<GzFastqReader> where P: AsRef<Path> {
//     let f = File::open(path)?;
//     Ok(FastqReader::new(MultiGzDecoder::new(BufReader::with_capacity(64 * 1024, f))))
// }

pub fn open<P>(path: P) -> io::Result<FastqReader<Box<dyn Read>>>
where
    P: AsRef<Path>,
{
    let f = File::open(&path)?;
    let buf = BufReader::with_capacity(64 * 1024, f);

    let reader: Box<dyn Read> = if path.as_ref().extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(buf))
    } else {
        Box::new(buf)
    };

    Ok(FastqReader::new(reader))
}