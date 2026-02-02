use crate::utils::position::Position;

use regex::Regex;

use clap::ValueEnum;

pub enum OpenstContext {
    Chip,
    Tile,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
pub enum BarcodeMode {
    Openst,
    SCMethyl,
    ChromiumATAC,
    ChromiumMRNA,
    Custom,
}

pub type BarcodeConfig = (Position, String);
impl BarcodeMode {
    pub fn openst(context: OpenstContext) -> BarcodeConfig {
        match context {
            OpenstContext::Chip => Self::openst_chip(),
            OpenstContext::Tile => Self::openst_tile(),
        }
    }

    fn openst_chip() -> BarcodeConfig {
        let pos = Position::new(true, false, true, 2, 30);
        // HDMI32-DraI: NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN
        let pattern = String::from("NNNBNNBNNBNNBNNBNNBNNBNNBVNB");
        (pos, pattern)
    }

    fn openst_tile() -> BarcodeConfig {
        let pos = Position::new(true, false, false, 2, 30);
        // barcode:     NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN
        let pattern = String::from("VNBVNNVNNVNNVNNVNNVNNVNNVNNN");
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
