#!/usr/bin/env python
"""Download gnomAD gene constraint scores (pLI, LOEUF)."""

from __future__ import annotations

import urllib.request
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "gnomad"

# gnomAD v4 constraint metrics
GNOMAD_CONSTRAINT_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/"
    "release/4.1/constraint/"
    "gnomad.v4.1.constraint_metrics.tsv"
)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    dest = OUTPUT_DIR / "gnomad.v4.1.constraint_metrics.tsv"

    if dest.exists():
        print(f"Already exists: {dest}")
        return

    print(f"Downloading gnomAD constraint metrics...")
    print(f"  URL: {GNOMAD_CONSTRAINT_URL}")
    print(f"  Destination: {dest}")

    try:
        urllib.request.urlretrieve(GNOMAD_CONSTRAINT_URL, dest)
        size_mb = dest.stat().st_size / 1024 / 1024
        print(f"  Done: {size_mb:.1f} MB")
    except Exception as e:
        print(f"  FAILED: {e}")
        # Try alternative URL (gnomAD v2.1.1 if v4 fails)
        alt_url = (
            "https://storage.googleapis.com/gcp-public-data--gnomad/"
            "release/2.1.1/constraint/"
            "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
        )
        alt_dest = OUTPUT_DIR / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
        print(f"  Trying v2.1.1 alternative: {alt_url}")
        try:
            urllib.request.urlretrieve(alt_url, alt_dest)
            print(f"  Done: {alt_dest.stat().st_size / 1024 / 1024:.1f} MB")
        except Exception as e2:
            print(f"  Also failed: {e2}")
            print("  Please download manually from https://gnomad.broadinstitute.org/downloads")


if __name__ == "__main__":
    main()
