
use crate::utils::{
    fastqfile::{open, FastqReader},
    position::Position,
    barcode_iter::{validate_absolute_dirpath, BarcodesIter},
    error::AppError,
};

use std::{fs, io::{self, BufWriter, Write}, process::Command};
use std::path::{PathBuf, Path};
use regex::Regex;
use clap::{Parser, ValueEnum};

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

#[derive(Parser, Debug)]
#[command(name = "bcl")]
#[command(about = "Process bcl dir into chip barcode list", long_about = None)]
#[command(next_line_help = true)]
pub struct TouchBarcodeArgs {
    /// Path to BCL directory
    #[arg(
        short = 'I', 
        long, 
        required = true,
        value_parser = validate_absolute_dirpath,
    )]
    bcl_dir: PathBuf,

    /// Path to output directory
    #[arg(short, long, required = true, value_parser = validate_absolute_dirpath)]
    output: PathBuf,

    /// barcode parsing mode
    #[arg(short, long, value_enum, default_value_t = BarcodeMode::Openst)]
    mode: BarcodeMode,

    /// turn on to run fastqc on each tile's fastq file
    #[arg(long)]
    fastqc: bool,

    /// Custom barcode position (only effective when mode=custom)
    /// 
    /// Format: "read{1/2}:{+/-}:start-end" 
    /// 
    /// Due to single-ended sequencing, there should only be read1, (e.g. "read1:+:1-16" or "read1:-:2-30")
    #[arg(
        long, 
        required_if_eq("mode", "custom"), 
        value_parser = clap::value_parser!(Position),
        value_name = "BARCODE_POS",
    )]
    barcode_pos: Option<Position>,

    /// Custom barcode pattern (only effective when mode=custom)
    /// 
    /// Regex: ^[ATGCNRYMKSWHBVD]+$
    /// 
    /// there should only be the pattern before convert sequence into reverse complement sequence.
    /// (e.g. openst-barcode: VNBVNNVNNVNNVNNVNNVNNVNNVNNN, openst-seq: NNNBNNBNNBNNBNNBNNBNNBNNBVNB)
    #[arg(
        long, 
        required_if_eq("mode", "custom"), 
        value_parser = validate_barcode_pattern,
        value_name = "BARCODE_PATTERN",
    )]
    barcode_pattern: Option<String>,
}

impl TouchBarcodeArgs {
    pub fn init(self) -> InitTouchBarcodeArgs {
        let (pos, pattern) = match (self.barcode_pos, self.barcode_pattern) {
            (Some(pos), Some(pattern)) => (pos, pattern),
            (None, None) => BarcodeMode::openst(),
            _ => unreachable!("clap parse the error is impossible.")
        };
        InitTouchBarcodeArgs::new(self.bcl_dir, self.output, self.fastqc, pos, pattern)
    }
}

pub struct InitTouchBarcodeArgs {
    bcl_dir: PathBuf,
    output: PathBuf,
    fastqc: bool,
    pos: Position,
    pattern: String,
}

impl InitTouchBarcodeArgs {
    #[inline]
    fn new(
        bcl_dir: PathBuf, 
        output: PathBuf, 
        fastqc: bool, 
        pos: Position, 
        pattern: String
    ) -> Self {
        Self {
            bcl_dir,
            output,
            fastqc,
            pos,
            pattern
        }
    }

    #[inline]
    fn bcl_dir(&self) -> &Path { self.bcl_dir.as_path() }

    #[inline]
    pub fn output(&self) -> &Path { &self.output.as_path() }

    #[inline]
    fn pos(&self) -> &Position { &self.pos }

    #[inline]
    fn pattern(&self) -> &str { &self.pattern }

    #[inline]
    pub fn fastq_path(&self, tile_id: &str) -> PathBuf { 
        self.output.join(format!("fastq/{tile_id}"))
    }

    #[inline]
    pub fn fastq_file(&self, tile_id: &str) -> PathBuf { 
        self.output.join(format!("fastq/{tile_id}/Undetermined_S0_R1_001.fastq.gz"))
    }

    #[inline]
    pub fn tmp_file(&self, tile_id: &str) -> PathBuf {
        self.output.join(format!("tmp/{}.txt", tile_id))
    }

    fn command_nonexists(&self, command: &str) -> io::Result<()> {
        let stauts = Command::new(command).arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .is_ok();
        if stauts {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("{} command not found", command),
            ))
        }
    }

    #[cfg(target_os = "macos")]
    fn docker_image_nonexists(&self, image: &str) -> io::Result<()> {
        let output = Command::new("docker").args(&["images", "-q", image]).output()?;

        if output.stdout.len() > 0 {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("{} image not found", image),
            ))
        }
    }

    pub fn validate_command(&self) -> io::Result<()> {
        if self.fastqc {
            self.command_nonexists("fastqc")?;
        }
        #[cfg(target_os = "linux")]
        self.command_nonexists("bcl-convert")?;
        #[cfg(target_os = "macos")]
        {
            self.command_nonexists("docker")?;
            self.docker_image_nonexists("zymoresearch/bcl-convert")?;
        }
        self.command_nonexists("bgzip")?;
        self.command_nonexists("tabix")
    }

    pub fn extract_tile_ids(&self) -> Result<Vec<String>, AppError> {
        let path = self.bcl_dir().join("RunInfo.xml");
        let re = Regex::new(r#"<Tile>([1-4]_[0-9]{4})</Tile>"#).unwrap();
        let content = fs::read_to_string(&path)?;
        let tile_ids: Vec<String> = re.captures_iter(&content)
        .filter_map(|cap| cap.get(1).map(
            |id| id.as_str().to_string()
        )).collect();
        if tile_ids.is_empty() { 
            return Err(AppError::EmptyTileIDsList(path)) 
        } else {
            Ok(tile_ids)
        }
    }

    fn run_command(
        &self,
        command: &str,
        args: &[&str],
        output_dir: &Path,
        tile_id: &str,
        error_msg: &str,
    ) -> Result<(), AppError> {
        use std::process::Stdio;
    
        // 确保输出目录存在
        if !output_dir.exists() {
            fs::create_dir_all(output_dir)?;
        }
        
        // 创建/打开日志文件（追加模式）
        let log_path = output_dir.join("command_output.log");
        let mut log_file = fs::OpenOptions::new().create(true).append(true).open(log_path)?;
        
        // 执行命令
        let output = Command::new(command).args(args)
            .stdout(Stdio::piped()).stderr(Stdio::piped()).output()?;
        
        // 记录日志
        writeln!(
            log_file,
            "{} stdout in tile_id {}:\n{}",
            command,
            tile_id,
            String::from_utf8_lossy(&output.stdout)
        )?;
        writeln!(
            log_file,
            "{} stderr in tile_id {}:\n{}",
            command,
            tile_id,
            String::from_utf8_lossy(&output.stderr)
        )?;
        
        // 检查执行状态
        if !output.status.success() {
            return Err(AppError::CommandError(
                format!("{} in tile_id {}", error_msg, tile_id)
            ));
        }
        
        Ok(())
    }

    fn bcl_convert(&self, tile_id: &str, fastq_dir: &Path) -> Result<(), AppError> {
        let args = [
            "--bcl-input-directory", &self.bcl_dir.display().to_string(),
            "--output-directory", &fastq_dir.display().to_string(),
            "--tiles", &format!("s_{}", tile_id),
            "--no-sample-sheet", "true",
            "--no-lane-splitting", "true",
            "--force"
        ];
        
        self.run_command(
            "bcl-convert",
            &args,
            &fastq_dir,
            tile_id,
            "bcl-convert run failed"
        )
    }
    
    fn docker_image_run(&self, tile_id: &str, fastq_dir: &Path) -> Result<(), AppError> {        
        let args = [
            "run", "--rm",
            "-v", &format!("{}:/mnt/run", self.bcl_dir.display()),
            "-v", &format!("{}:/mnt/output", fastq_dir.display()),
            "zymoresearch/bcl-convert",
            "--bcl-input-directory", "/mnt/run",
            "--output-directory", "/mnt/output",
            "--tiles", &format!("s_{}", tile_id),
            "--no-sample-sheet", "true",
            "--no-lane-splitting", "true",
            "--force"
        ];
        
        self.run_command(
            "docker",
            &args,
            &fastq_dir,
            tile_id,
            "Docker run failed"
        )
    }

    fn fastqc_run(&self, tile_id: &str) -> Result<(), AppError> {
        let fastq_file = self.fastq_file(tile_id);
        
        self.run_command(
            "fastqc",
            &[fastq_file.as_os_str().to_str().unwrap()],
            &self.fastq_path(tile_id),
            tile_id,
            "FastQC failed"
        )
    }

    pub fn convert_bcl_into_tile(&self, tile_id: &str) -> Result<(), AppError> {
        let fastq_dir = self.fastq_path(tile_id);
        if cfg!(target_os = "linux") {
            self.bcl_convert(tile_id, &fastq_dir)?;
        } else if cfg!(target_os = "macos") {
            self.docker_image_run(tile_id, &fastq_dir)?;
        } else {
            return Err(AppError::UnsupportedOS);
        }
    
        if self.fastqc {
            self.fastqc_run(tile_id)?;
        }
        Ok(())
    }

    pub fn create_barcode_iter(&self, tile_id: &str) -> io::Result<BarcodesIter<BufWriter<fs::File>>> {
        let inner: FastqReader = open(
            self.fastq_path(tile_id).join("Undetermined_S0_R1_001.fastq.gz")
        )?;
        let tmp_path = self.tmp_file(tile_id);
        let writer = fs::OpenOptions::new().write(true)
            .create(true).open(tmp_path).map(BufWriter::new)?;
        Ok(BarcodesIter::into_file(inner, self.pos(), self.pattern(), writer))
    }
}


#[derive(ValueEnum, Clone, Copy, Debug)]
enum BarcodeMode {
    Openst,
    Custom,
}

pub type BarcodeConfig = (Position, String);
impl BarcodeMode {
    pub fn openst() -> BarcodeConfig {
        let pos = Position::new(false, true, 2, 30);
        // HDMI32-DraI: NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN
        // revcomp:     NNNNNBNNBNNBNNBNNBNNBNNBNNBVNBNN
        let pattern: String = String::from("NNNBNNBNNBNNBNNBNNBNNBNNBVNB");
        (pos, pattern)
    }
}