"""Download AlphaFold PPI predictions from Dryad (Baker Lab, RoseTTAFold2-PPI).

Source: Zhang et al., Science 2024 (doi:10.1126/science.adt1630)
Data: https://datadryad.org/dataset/doi:10.5061/dryad.15dv41p84

Downloads final_predictions.tar.gz containing:
  - final_predictions_90.tsv (17,849 interactions at 90% precision)
  - final_predictions_80.tsv (broader set at 80% precision)
"""

from __future__ import annotations

import tarfile
import urllib.request
from pathlib import Path

# Dryad file stream ID for final_predictions.tar.gz (14.13 MB)
DRYAD_FILE_URL = "https://datadryad.org/downloads/file_stream/4357020"

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "raw" / "networks" / "alphafold_baker"


def download_alphafold_ppi(output_dir: Path | None = None) -> Path:
    """Download and extract AlphaFold PPI predictions.

    Returns:
        Path to the extracted directory containing TSV files.
    """
    out = output_dir or DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    target_file = out / "final_predictions_90.tsv"
    if target_file.exists():
        print(f"Already exists: {target_file}")
        return target_file

    archive = out / "final_predictions.tar.gz"
    print(f"Downloading AlphaFold PPI predictions from Dryad...")
    print(f"  URL: {DRYAD_FILE_URL}")
    print(f"  Size: ~14 MB")
    req = urllib.request.Request(
        DRYAD_FILE_URL,
        headers={"User-Agent": "SSA-ME/2.0 (research; https://github.com)"},
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        archive.write_bytes(resp.read())
    print(f"  Saved: {archive}")

    print("Extracting...")
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(path=out, filter="data")
    print(f"  Extracted to: {out}")

    # Find the extracted TSV
    for tsv in out.rglob("final_predictions_90.tsv"):
        print(f"  Found: {tsv}")
        return tsv

    # Fallback: list what was extracted
    extracted = list(out.rglob("*.tsv"))
    if extracted:
        print(f"  Found TSV files: {[str(f.name) for f in extracted]}")
        return extracted[0]

    print(f"  Warning: No TSV files found in extracted archive.")
    return out


if __name__ == "__main__":
    result = download_alphafold_ppi()
    print(f"\nDone. Data at: {result}")
