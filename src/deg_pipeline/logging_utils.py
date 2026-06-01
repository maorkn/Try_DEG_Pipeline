"""Logging and manifest utilities for the DEG pipeline."""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import sys
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


def configure_logging(outdir: Path, script_name: str | None = None) -> logging.Logger:
    """Configure logging for a pipeline script.
    
    Args:
        outdir: Output directory where logs will be written.
        script_name: Optional script name for log file naming.
        
    Returns:
        Configured logger instance.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    
    if script_name is None:
        script_name = Path(sys.argv[0]).stem
    
    log_file = outdir / f"{script_name}.log"
    
    # Configure root logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    
    logger = logging.getLogger(script_name)
    logger.info(f"Logging to {log_file}")
    return logger


def compute_file_checksum(filepath: Path, algorithm: str = "sha256") -> str:
    """Compute checksum of a file.
    
    Args:
        filepath: Path to the file.
        algorithm: Hash algorithm to use (default: sha256).
        
    Returns:
        Hexadecimal checksum string.
    """
    if filepath.is_dir():
        return "directory"

    hasher = hashlib.new(algorithm)
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def get_git_commit_hash() -> str | None:
    """Get the current git commit hash if available.
    
    Returns:
        Git commit hash string or None if not in a git repository.
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_python_packages() -> dict[str, str]:
    """Get versions of installed Python packages relevant to the pipeline.
    
    Returns:
        Dictionary mapping package names to version strings.
    """
    relevant_packages = [
        "pandas",
        "numpy",
        "scipy",
        "statsmodels",
        "pydeseq2",
        "pyyaml",
        "pydantic",
        "matplotlib",
        "seaborn",
        "scikit-learn",
        "jinja2",
        "goatools",
        "deg-pipeline",
    ]
    
    versions = {}
    for pkg_name in relevant_packages:
        try:
            versions[pkg_name] = version(pkg_name)
        except PackageNotFoundError:
            versions[pkg_name] = "not installed"
    
    return versions


def write_manifest(
    outdir: Path,
    script_name: str,
    config: dict[str, Any] | None = None,
    input_files: list[Path] | None = None,
    output_files: list[Path] | None = None,
    extra_info: dict[str, Any] | None = None,
    filename: str = "manifest.json",
) -> Path:
    """Write a reproducibility manifest for a pipeline stage.
    
    Args:
        outdir: Output directory for the manifest.
        script_name: Name of the script that produced these outputs.
        config: Configuration dictionary used for this run.
        input_files: List of input file paths to checksum.
        output_files: List of output file paths to checksum.
        extra_info: Additional information to include in the manifest.
        filename: Manifest filename (default: manifest.json).
        
    Returns:
        Path to the written manifest file.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / filename
    
    manifest: dict[str, Any] = {
        "script": script_name,
        "timestamp_start": datetime.now().isoformat(),
        "python_version": sys.version,
        "git_commit": get_git_commit_hash(),
        "packages": get_python_packages(),
    }
    
    if config is not None:
        manifest["config"] = config
    
    if input_files:
        manifest["input_files"] = {
            str(f): compute_file_checksum(f) if f.exists() else "file not found"
            for f in input_files
        }
    
    if output_files:
        manifest["output_files"] = {
            str(f): compute_file_checksum(f) if f.exists() else "file not found"
            for f in output_files
        }
    
    if extra_info:
        manifest.update(extra_info)
    
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    
    return manifest_path


def check_output_exists(outdir: Path, filename: str, force: bool = False) -> Path | None:
    """Check if an output file already exists.
    
    Args:
        outdir: Output directory.
        filename: Output filename.
        force: If True, return None to allow overwriting.
        
    Returns:
        Path to existing file if it exists and force is False, else None.
        
    Raises:
        FileExistsError: If file exists and force is False.
    """
    filepath = outdir / filename
    if filepath.exists() and not force:
        raise FileExistsError(
            f"Output file {filepath} already exists. Use --force to overwrite."
        )
    return None if force or not filepath.exists() else filepath
