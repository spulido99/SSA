#!/usr/bin/env python
"""Download DepMap CRISPR essentiality data."""

from __future__ import annotations

import urllib.request
from pathlib import Path

# DepMap Public 24Q2 release (via Figshare - stable download URLs)
# Article: https://plus.figshare.com/articles/dataset/DepMap_24Q2_Public/25880521
FILES = {
    "CRISPRGeneEffect.csv": "https://ndownloader.figshare.com/files/46489063",  # ~400 MB
    "Model.csv": "https://ndownloader.figshare.com/files/46489732",  # ~545 KB
}

OUTPUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "depmap"


def download_file(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"  Already exists: {dest.name}")
        return
    print(f"  Downloading {dest.name}...")
    urllib.request.urlretrieve(url, dest)
    print(f"  Done: {dest.stat().st_size / 1024 / 1024:.1f} MB")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")

    for filename, url in FILES.items():
        download_file(url, OUTPUT_DIR / filename)

    print("\nDepMap download complete.")
    print("Files:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        print(f"  {f.name}: {f.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
