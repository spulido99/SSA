"""Download IntOGen driver gene catalog.

Source: IntOGen (intogen.org), latest release 2024.09.20
Data: ~600 driver genes from 28K tumors, 66 cancer types
License: Open access
"""

from __future__ import annotations

import zipfile
from pathlib import Path
from urllib.request import urlretrieve

INTOGEN_URL = "https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "raw" / "ground_truth"


def download_intogen(output_dir: Path | None = None) -> Path:
    """Download and extract IntOGen drivers catalog.

    Returns:
        Path to drivers.tsv file.
    """
    out = output_dir or DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    target = out / "intogen" / "drivers.tsv"
    if target.exists():
        print(f"Already exists: {target}")
        return target

    archive = out / "IntOGen-Drivers-20240920.zip"
    print("Downloading IntOGen driver catalog...")
    print(f"  URL: {INTOGEN_URL}")
    urlretrieve(INTOGEN_URL, archive)
    print(f"  Saved: {archive}")

    extract_dir = out / "intogen"
    extract_dir.mkdir(exist_ok=True)
    print("Extracting...")
    with zipfile.ZipFile(archive, "r") as zf:
        zf.extractall(extract_dir)
    print(f"  Extracted to: {extract_dir}")

    # Prefer Compendium_Cancer_Genes.tsv (canonical filtered drivers)
    # over Unfiltered_drivers.tsv (includes non-driver entries)
    for preferred in ("Compendium_Cancer_Genes.tsv", "Unfiltered_drivers.tsv"):
        for tsv in extract_dir.rglob(preferred):
            print(f"  Found: {tsv}")
            return tsv

    # Fallback: any TSV with "driver" in the name
    for tsv in extract_dir.rglob("*.tsv"):
        if "driver" in tsv.name.lower():
            print(f"  Found: {tsv}")
            return tsv

    return extract_dir


if __name__ == "__main__":
    result = download_intogen()
    print(f"\nDone. Data at: {result}")
