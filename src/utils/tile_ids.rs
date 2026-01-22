pub fn is_valid_tile_id(value: &str) -> Result<u64, String> {
    let tile_id: u64 = value.parse().map_err(|_| format!("`{}` is not valid integer", value))?;

    if VALID_TILE_IDS.contains(&tile_id) {
        Ok(tile_id)
    } else {
        Err(format!("tile_id {} is not in the valid range (valid range: 11101-42678)", tile_id))
    }
}

pub const VALID_TILE_IDS: [u64; 3744] = {
    // Array size: 4 × 2 × 6 × 78 = 3744
    let mut result = [0u64; 3744];
    let mut index = 0;

    let mut a = 1;
    while a <= 4 {
        let mut b = 1;
        while b <= 2 {
            let mut c = 1;
            while c <= 6 {
                let mut d = 1;
                while d <= 78 {
                    result[index] = a * 10000 + b * 1000 + c * 100 + d;
                    index += 1;
                    d += 1;
                }
                c += 1;
            }
            b += 1;
        }
        a += 1;
    }
    result
};