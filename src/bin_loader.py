import os
import sys
import platform
import subprocess
from pathlib import Path


def get_platform():
    """Returns the platform"""
    system = platform.system()
    arch = platform.machine().lower()

    if system == "Darwin":
        if arch == "x86_64":
            return "x86_64-apple-darwin"
        elif arch in ("arm64", "aarch64"):
            return "aarch64-apple-darwin"

    elif system == "Linux":
        try:
            output = subprocess.check_output(
                ["ldd", "--version"], stderr=subprocess.STDOUT
            )
            if b"musl" in output:
                return "x86_64-unknown-linux-musl"
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass

        if arch == "x86_64":
            return "x86_64-unknown-linux-gnu"
        
    return None

def get_binary_path():
    """Returns the binary file path depending on the platform"""
    platform_id = get_platform()
    if not platform_id:
        raise RuntimeError(
            f"Unsupported platform: {platform.system()}/{platform.machine()}"
        )

    bin_path = Path(__file__).parent / "binaries" / platform_id / "opentools"

    if not bin_path.exists():
        dev_path = Path.cwd() / "binaries" / platform_id / "opentools"
        if dev_path.exists():
            return dev_path

        raise FileNotFoundError(f"Binary not found for platform {platform_id}")

    return bin_path

def run_binary():
    """Execute the binary"""
    bin_path = get_binary_path()
    bin_path = bin_path.resolve()

    try:
        os.chmod(bin_path, 0o755)
    except Exception as e:
        sys.stderr.write(f"Warning: Could not set executable permissions: {e}\n")

    result = subprocess.run([str(bin_path)] + sys.argv[1:])
    sys.exit(result.returncode)
