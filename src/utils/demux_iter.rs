use std::io::{ self, BufWriter, Read };
use std::fs::File;

use seq_io::fastq::{ OwnedRecord, Record, Reader as FastqReader };

use crate::utils::{
    barcode_iter::BarcodesIter,
    error::AppError,
    position::Position
};

pub struct DemuxIter<'a> {
    r1: FastqReader<Box<dyn Read>>,
    r2: FastqReader<Box<dyn Read>>,
    i2: BarcodesIter<'a, Box<dyn Read>>,
}

impl<'a> DemuxIter<'a> {
    #[inline]
    pub fn new(r1: FastqReader<Box<dyn Read>>, r2: FastqReader<Box<dyn Read>>, i2: BarcodesIter<'a, Box<dyn Read>>) -> Self {
        Self { r1, r2, i2 }
    }

    #[inline]
    pub fn pos(&self) -> &Position {
        self.i2.pos()
    }

    pub fn for_each<F>(&mut self, mut f: F) -> Result<(), AppError>
        where F: FnMut(DemuxRead) -> Result<(), AppError>
    {
        loop {
            let r1 = match self.r1.records().next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => {
                    return Err(e.into());
                }
                None => {
                    break;
                }
            };

            let r2 = match self.r2.records().next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => {
                    return Err(e.into());
                }
                None => {
                    break;
                }
            };

            let bc = match self.i2.next() {
                Some(Ok(bc)) => bc,
                Some(Err(e)) => {
                    return Err(e);
                }
                None => {
                    break;
                }
            };

            f(DemuxRead::new(r1, r2, bc))?;
        }
        Ok(())
    }
}

impl<'a> Iterator for DemuxIter<'a> {
    type Item = Result<DemuxRead, AppError>;

    fn next(&mut self) -> Option<Self::Item> {
        let r1 = match self.r1.records().next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => {
                return Some(Err(e.into()));
            }
            None => {
                return None;
            }
        };

        let r2 = match self.r2.records().next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => {
                return Some(Err(e.into()));
            }
            None => {
                return None;
            }
        };

        let i2 = match self.i2.next() {
            Some(Ok(i2)) => i2,
            Some(Err(e)) => {
                return Some(Err(e));
            }
            None => {
                return None;
            }
        };

        Some(Ok(DemuxRead::new(r1, r2, i2)))
    }
}

pub struct DemuxRead {
    r1: OwnedRecord,
    r2: OwnedRecord,
    bc: String,
}

impl DemuxRead {
    #[inline]
    pub fn new(r1: OwnedRecord, r2: OwnedRecord, bc: String) -> Self {
        Self {
            r1,
            r2,
            bc,
        }
    }

    #[inline]
    pub fn read1_head(&self) -> &[u8] {
        self.r1.head()
    }

    #[inline]
    pub fn read1_seq(&self) -> &[u8] {
        self.r1.seq()
    }

    #[inline]
    pub fn read1_qual(&self) -> &[u8] {
        self.r1.qual()
    }

    #[inline]
    pub fn read2_head(&self) -> &[u8] {
        self.r2.head()
    }

    #[inline]
    pub fn read2_seq(&self) -> &[u8] {
        self.r2.seq()
    }

    #[inline]
    pub fn read2_qual(&self) -> &[u8] {
        self.r2.qual()
    }

    #[inline]
    pub fn barcode(&self) -> &str {
        &self.bc
    }

    pub fn write_r1<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        self.write_with_bc(&self.r1, writer)
    }

    pub fn write_r2<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        self.write_with_bc(&self.r1, writer)
    }

    pub fn write_i2<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        let mut rec = self.r1.clone();
        let qual = vec![b'I'; self.barcode().len()];

        rec.seq = self.barcode().to_string().into_bytes();
        rec.qual = qual;

        self.write_with_bc(&rec, writer)
    }

    fn write_with_bc<W: io::Write>(&self, rec: &OwnedRecord, writer: &mut W) -> io::Result<()> {
        writer.write_all(b"@")?;
        writer.write_all(rec.head())?;
        writer.write_all(b"\n")?;
        writer.write_all(rec.seq())?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(rec.qual())?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}

pub struct DemuxWriters {
    pub r1: BufWriter<File>,
    pub r2: BufWriter<File>,
    pub i2: BufWriter<File>,
}
