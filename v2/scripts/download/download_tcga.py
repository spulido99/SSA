#!/usr/bin/env python
"""Download TCGA MC3 somatic mutation data."""

from __future__ import annotations

import urllib.request
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "tcga"

# MC3 MAF: controlled-access via GDC, but public summary available
# The full MC3 MAF is at GDC but requires dbGaP authorization.
# We use the open-access mutation counts and TCGA PanCanAtlas summary.

DOWNLOADS = {
    # TCGA PanCanAtlas: MC3 public MAF (open access subset)
    "mc3.v0.2.8.PUBLIC.maf.gz": (
        "https://api.gdc.cancer.gov/data/"
        "1c8cfe5f-e52d-41ba-94da-f15ea1337efc"
    ),
}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")

    for filename, url in DOWNLOADS.items():
        dest = OUTPUT_DIR / filename
        if dest.exists():
            print(f"Already exists: {dest.name}")
            continue

        print(f"Downloading {filename}...")
        print(f"  URL: {url}")
        print(f"  This is ~50MB, may take a few minutes...")

        try:
            urllib.request.urlretrieve(url, dest)
            size_mb = dest.stat().st_size / 1024 / 1024
            print(f"  Done: {size_mb:.1f} MB")
        except Exception as e:
            print(f"  FAILED: {e}")
            print("  The MC3 MAF requires GDC access.")
            print("  Alternative: download from https://gdc.cancer.gov/about-data/publications/mc3-2017")
            if dest.exists():
                dest.unlink()

    print("\nTCGA download complete.")


if __name__ == "__main__":
    main()
