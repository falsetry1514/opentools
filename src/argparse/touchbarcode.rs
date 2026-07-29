// ------------------
// 1. std
// ------------------
use std::fs::{ self, File };
use std::io::{ self, BufWriter, Write };
use std::path::{ Path, PathBuf };
use std::process::Command;

// ------------------
// 2. external crates
// ------------------
use clap::Parser;
use regex::Regex;
use sysinfo::System;

// ------------------
// 3. inner crate
// ------------------
use crate::utils::{
    barcode_config::{ BarcodeMode, OpenstContext, validate_barcode_pattern },
    barcode_iter::{ BarcodesIter, validate_absolute_dirpath },
    error::AppError,
    fastq_iter::open,
    position::Position,
};

#[derive(Parser, Debug)]
#[command(name = "bcl")]
#[command(about = "Process bcl dir into chip barcode list", long_about = None)]
#[command(next_line_help = true)]
pub struct TouchBarcodeArgs {
    /// Path to BCL directory
    #[arg(value_parser = validate_absolute_dirpath)]
    bcl_dir: PathBuf,

    /// Path to output directory
    /// If not provided, defaults to the parent of `bcl_dir`.
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,

    /// The cores used for running
    #[arg(short = '@', long)]
    cores: Option<usize>,

    /// barcode parsing mode
    #[arg(short, long, value_enum, default_value_t = BarcodeMode::Openst28bc)]
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
        value_name = "BARCODE_POS"
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
        value_name = "BARCODE_PATTERN"
    )]
    barcode_pattern: Option<String>,
}

impl TouchBarcodeArgs {
    pub fn init(self) -> Result<InitTouchBarcodeArgs, AppError> {
        let output = self.output.unwrap_or_else(|| {
            self.bcl_dir.parent().expect("bcl_dir has no parent").to_path_buf()
        });

        let (pos, pattern) = match (self.mode, self.barcode_pos, self.barcode_pattern) {
            (BarcodeMode::Custom, Some(pos), Some(pattern)) => (pos, pattern),
            (BarcodeMode::Custom, _, _) => {
                return Err(
                    AppError::CommandError(
                        "Custom barcode mode requires --barcode-pos and --barcode-pattern".into()
                    )
                );
            }
            (BarcodeMode::Openst28bc, None, None) => BarcodeMode::openst(OpenstContext::Chip28),
            (BarcodeMode::Openst32bc, None, None) => BarcodeMode::openst(OpenstContext::Chip32),
            _ => {
                unimplemented!("Other barcode modes are unimplemented!");
            }
        };

        let upper_cores = {
            let mut sys = System::new();
            sys.refresh_cpu_all();
            sys.cpus().len().saturating_sub(2).max(1)
        };

        let n_core = self.cores.unwrap_or(upper_cores).min(upper_cores);

        Ok(InitTouchBarcodeArgs::new(self.bcl_dir, output, n_core, self.fastqc, pos, pattern))
    }
}

pub struct InitTouchBarcodeArgs {
    bcl_dir: PathBuf,
    output: PathBuf,
    cores: usize,
    fastqc: bool,
    pos: Position,
    pattern: String,
}

impl InitTouchBarcodeArgs {
    #[inline]
    fn new(
        bcl_dir: PathBuf,
        output: PathBuf,
        cores: usize,
        fastqc: bool,
        pos: Position,
        pattern: String
    ) -> Self {
        Self {
            bcl_dir,
            output,
            cores,
            fastqc,
            pos,
            pattern,
        }
    }

    #[inline]
    fn bcl_dir(&self) -> &Path {
        self.bcl_dir.as_path()
    }

    #[inline]
    pub fn output(&self) -> &Path {
        &self.output.as_path()
    }

    #[inline]
    pub fn cores(&self) -> usize {
        self.cores
    }

    #[inline]
    fn pos(&self) -> &Position {
        &self.pos
    }

    #[inline]
    fn pattern(&self) -> &str {
        &self.pattern
    }

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
        Command::new(command)
            .arg("--version")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map(|_| ()) // 如果成功，返回 ()
            .map_err(|_|
                io::Error::new(io::ErrorKind::NotFound, format!("{} command not found", command))
            )
    }

    #[cfg(target_os = "macos")]
    fn docker_image_nonexists(&self, image: &str) -> io::Result<()> {
        Command::new("docker")
            .args(&["images", "-q", image])
            .output()
            .map(|output| {
                if !output.stdout.is_empty() {
                    Ok(())
                } else {
                    Err(io::Error::new(io::ErrorKind::NotFound, format!("{image} image not found")))
                }
            })?
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
        let tile_ids: Vec<String> = re
            .captures_iter(&content)
            .filter_map(|cap| cap.get(1).map(|id| id.as_str().to_string()))
            .collect();
        (!tile_ids.is_empty()).then(|| tile_ids).ok_or_else(|| AppError::EmptyTileIDsList(path))
    }

    fn run_command(
        &self,
        command: &str,
        args: &[&str],
        output_dir: &Path,
        tile_id: &str,
        error_msg: &str
    ) -> Result<(), AppError> {
        use std::process::Stdio;

        fs::create_dir_all(output_dir)?;

        // 创建/打开日志文件（追加模式）
        let log_path = output_dir.join("command_output.log");
        let mut log_file = fs::OpenOptions
            ::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(log_path)?;

        // 执行命令
        let output = Command::new(command)
            .args(args)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()?;

        // 记录日志
        for (stream_name, buf) in &[
            ("stdout", &output.stdout),
            ("stderr", &output.stderr),
        ] {
            writeln!(
                log_file,
                "{} {} in tile_id {}:\n{}",
                command,
                stream_name,
                tile_id,
                String::from_utf8_lossy(buf)
            )?;
        }

        // 检查执行状态
        output.status.success()
            .then(|| ())
            .ok_or_else(|| AppError::CommandError(format!("{} in tile_id {}", error_msg, tile_id)))
    }

    fn bcl_convert(&self, tile_id: &str, fastq_dir: &Path) -> Result<(), AppError> {
        let args = [
            "--bcl-input-directory",
            &self.bcl_dir.display().to_string(),
            "--output-directory",
            &fastq_dir.display().to_string(),
            "--tiles",
            &format!("s_{}", tile_id),
            "--no-sample-sheet",
            "true",
            "--no-lane-splitting",
            "true",
            "--force",
        ];

        self.run_command("bcl-convert", &args, &fastq_dir, tile_id, "bcl-convert run failed")
    }

    fn docker_image_run(&self, tile_id: &str, fastq_dir: &Path) -> Result<(), AppError> {
        let args = [
            "run",
            "--rm",
            "-v",
            &format!("{}:/mnt/run", self.bcl_dir.display()),
            "-v",
            &format!("{}:/mnt/output", fastq_dir.display()),
            "zymoresearch/bcl-convert",
            "--bcl-input-directory",
            "/mnt/run",
            "--output-directory",
            "/mnt/output",
            "--tiles",
            &format!("s_{}", tile_id),
            "--no-sample-sheet",
            "true",
            "--no-lane-splitting",
            "true",
            "--force",
        ];

        self.run_command("docker", &args, &fastq_dir, tile_id, "Docker run failed")
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
        match () {
            _ if cfg!(target_os = "linux") => self.bcl_convert(tile_id, &fastq_dir)?,
            _ if cfg!(target_os = "macos") => self.docker_image_run(tile_id, &fastq_dir)?,
            _ => return Err(AppError::UnsupportedOS),
        }

        if self.fastqc {
            self.fastqc_run(tile_id)?;
        }
        Ok(())
    }

    pub fn render_writer(&self, tile_id: &str) -> io::Result<BufWriter<File>> {
        fs::OpenOptions
            ::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(self.tmp_file(tile_id))
            .map(BufWriter::new)
    }

    pub fn create_barcode_iter(
        &self,
        tile_id: &str
    ) -> io::Result<BarcodesIter<'_, impl std::io::Read>> {
        let inner = open(self.fastq_path(tile_id).join("Undetermined_S0_R1_001.fastq.gz"))?;
        Ok(BarcodesIter::new(inner, self.pos(), self.pattern()))
    }
}

pub struct TouchBarcodeReport {
    total_count: u64,
    filter_qual_count: u64,
    filter_seq_count: u64,
    filter_dup_count: u64,
}

impl TouchBarcodeReport {
    #[inline]
    pub fn new(
        total_count: u64,
        filter_qual_count: u64,
        filter_seq_count: u64,
        filter_dup_count: u64
    ) -> Self {
        Self {
            total_count,
            filter_qual_count,
            filter_seq_count,
            filter_dup_count,
        }
    }

    #[inline]
    fn filtered_count(&self) -> u64 {
        self.filter_qual_count + self.filter_seq_count + self.filter_dup_count
    }

    #[inline]
    fn passed_count(&self) -> u64 {
        self.total_count - self.filtered_count()
    }
}

impl std::fmt::Display for TouchBarcodeReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Total={}, Filtered={} (Qual={}, Seq={}, Dup={}), Passed={}",
            self.total_count,
            self.filtered_count(),
            self.filter_qual_count,
            self.filter_seq_count,
            self.filter_dup_count,
            self.passed_count()
        )
    }
}
