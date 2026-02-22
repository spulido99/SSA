#!/usr/bin/env python
"""Download PPI network files: HI-II-14, Reactome, STRING."""

from __future__ import annotations

import gzip
import shutil
import ssl
import urllib.request
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "networks"

DOWNLOADS = {
    # HI-II-14: Human Interactome (unbiased Y2H screen)
    "HI-II-14.tsv": "http://www.interactome-atlas.org/data/HI-II-14.tsv",
    # HuRI: Human Reference Interactome (updated, larger unbiased screen)
    "HuRI.tsv": "http://www.interactome-atlas.org/data/HuRI.tsv",
    # Reactome: Functional interactions (GMT format)
    "ReactomePathways.gmt": (
        "https://reactome.org/download/current/ReactomePathways.gmt.zip"
    ),
    # STRING: protein links (combined score, human only, gzipped)
    "9606.protein.links.v12.0.txt.gz": (
        "https://stringdb-downloads.org/download/"
        "protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
    ),
    # STRING: aliases for ENSP -> gene symbol mapping
    "9606.protein.aliases.v12.0.txt.gz": (
        "https://stringdb-downloads.org/download/"
        "protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
    ),
}


def _ssl_context() -> ssl.SSLContext:
    """SSL context that skips certificate verification (interactome-atlas cert expired)."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def download_file(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"  Already exists: {dest.name}")
        return
    print(f"  Downloading {dest.name}...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "SSA-ME/2.0"})
        with urllib.request.urlopen(req, context=_ssl_context(), timeout=120) as resp:
            dest.write_bytes(resp.read())
        size_mb = dest.stat().st_size / 1024 / 1024
        print(f"  Done: {size_mb:.1f} MB")
    except Exception as e:
        print(f"  FAILED: {e}")
        if dest.exists():
            dest.unlink()


def decompress_gz(gz_path: Path) -> None:
    """Decompress .gz file, keeping original."""
    out_path = gz_path.with_suffix("")
    if out_path.exists():
        print(f"  Already decompressed: {out_path.name}")
        return
    print(f"  Decompressing {gz_path.name}...")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"  Done: {out_path.stat().st_size / 1024 / 1024:.1f} MB")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")

    for filename, url in DOWNLOADS.items():
        download_file(url, OUTPUT_DIR / filename)

    # Decompress STRING files
    for gz_file in OUTPUT_DIR.glob("*.gz"):
        decompress_gz(gz_file)

    print("\nNetwork download complete.")
    print("Files:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        if not f.name.endswith(".gz"):
            print(f"  {f.name}: {f.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
