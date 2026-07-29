use crate::utils::position::Position;

use regex::Regex;

use clap::ValueEnum;

pub enum OpenstContext {
    Chip28,
    Tile28,
    Chip32,
    Tile32,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
pub enum BarcodeMode {
    Openst28bc,
    Openst32bc,
    SCMethyl,
    ChromiumATAC,
    ChromiumMRNA,
    Custom,
}
const OPENST_CHIP_PATTERN: &[u8] = b"NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN";
const OPENST_TILE_PATTERN: &[u8] = b"NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN";

pub type BarcodeConfig = (Position, String);
impl BarcodeMode {
    pub fn openst(context: OpenstContext) -> BarcodeConfig {
        match context {
            OpenstContext::Chip28 => Self::openst_chip(2, 30),
            OpenstContext::Tile28 => Self::openst_tile(2, 30),
            OpenstContext::Chip32 => Self::openst_chip(0, 32),
            OpenstContext::Tile32 => Self::openst_tile(0, 32),
        }
    }

    fn openst_chip(start: usize, end: usize) -> BarcodeConfig {
        let pos = Position::new(true, false, true, start, end);
        // HDMI32-DraI: NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN
        let pattern = String::from_utf8_lossy(&OPENST_CHIP_PATTERN[start..end]).to_string();
        (pos, pattern)
    }

    fn openst_tile(start: usize, end: usize) -> BarcodeConfig {
        let pos = Position::new(true, false, false, start, end);
        // barcode:     NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN
        let pattern = String::from_utf8_lossy(&OPENST_TILE_PATTERN[start..end]).to_string();
        (pos, pattern)
    }
    
    pub fn chromium_mrna() -> BarcodeConfig {
        let pos = Position::new(true, false, false, 0, 16);
        let pattern = String::from("NNNNNNNNNNNNNNNN");
        (pos, pattern)
    }
    
    pub fn chromium_atac() -> BarcodeConfig {
        let pos = Position::new(false, true, true, 8, 24);
        let pattern = String::from("NNNNNNNNNNNNNNNN");
        (pos, pattern)
    }

    pub fn sc_methyl() -> BarcodeConfig {
        let pos = Position::new(false, true, true, 8, 24);
        let pattern = String::from("HHHHHHHHHHHHHHHH");
        (pos, pattern)
    }

}

pub fn validate_barcode_pattern(s: &str) -> Result<String, String> {
    let re = Regex::new(r"^[ATGCURYMKSWHBVDN]+$").unwrap();
    if re.is_match(s) {
        Ok(s.to_string())
    } else {
        Err(
            "Invalid barcode pattern. 
            Allowed characters: A, T, G, C, R, Y, M, K, S, W, H, B, V, D, N".to_string()
        )
    }
}
