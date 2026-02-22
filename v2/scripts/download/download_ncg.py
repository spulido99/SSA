"""Download NCG (Network of Cancer Genes) driver gene list.

Source: NCG 7.1 (Genome Biology 2022, doi:10.1186/s13059-022-02607-z)
Data: 3,177 cancer drivers with annotations
License: Open access (CC BY 4.0)

Strategy:
1. Try NCG website directly (unreliable — SSL cert expired / download.php broken)
2. Fall back to Springer supplementary data (Table S5 + Table S3 from PMC)
"""

from __future__ import annotations

import ssl
import urllib.request
from pathlib import Path

import pandas as pd

# NCG website URLs (often broken)
NCG_URLS = [
    "https://ncg.kcl.ac.uk/download.php?file=cancerdrivers_annotation",
    "http://network-cancer-genes.org/download.php?file=cancerdrivers_annotation",
]

# Springer supplementary data (stable)
SPRINGER_BASE = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1186%2Fs13059-022-02607-z/MediaObjects"
)
TABLE_S5_URL = f"{SPRINGER_BASE}/13059_2022_2607_MOESM6_ESM.xlsx"  # All drivers
TABLE_S3_URL = f"{SPRINGER_BASE}/13059_2022_2607_MOESM4_ESM.xlsx"  # Canonical + roles

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "raw" / "ground_truth"


def _try_ncg_website(target: Path) -> bool:
    """Try downloading from NCG website. Returns True if successful."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    for url in NCG_URLS:
        try:
            print(f"  Trying NCG website: {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "SSA-ME/2.0"})
            with urllib.request.urlopen(req, context=ctx, timeout=30) as resp:
                data = resp.read()
                # Verify it's actually TSV, not HTML
                if data[:10].startswith(b"<") or b"<!DOCTYPE" in data[:100]:
                    print("    Got HTML instead of TSV — website is broken")
                    continue
                if len(data) > 1000:
                    target.write_bytes(data)
                    print(f"    Saved: {target} ({len(data)} bytes)")
                    return True
                print(f"    Response too small ({len(data)} bytes)")
        except Exception as e:
            print(f"    Failed: {e}")
    return False


def _build_from_springer(output_dir: Path) -> Path:
    """Build NCG TSV from Springer supplementary data."""
    print("  Downloading supplementary data from Springer...")

    headers = {"User-Agent": "Mozilla/5.0 (compatible; SSA-ME/2.0)"}

    # Download Table S5 (all drivers)
    s5_path = output_dir / "NCG_TableS5_all_drivers.xlsx"
    if not s5_path.exists():
        req = urllib.request.Request(TABLE_S5_URL, headers=headers)
        with urllib.request.urlopen(req, timeout=60) as resp:
            s5_path.write_bytes(resp.read())
        print(f"    Table S5: {s5_path.stat().st_size:,} bytes")

    # Download Table S3 (canonical + roles)
    s3_path = output_dir / "NCG_TableS3_canonical_drivers.xlsx"
    if not s3_path.exists():
        req = urllib.request.Request(TABLE_S3_URL, headers=headers)
        with urllib.request.urlopen(req, timeout=60) as resp:
            s3_path.write_bytes(resp.read())
        print(f"    Table S3: {s3_path.stat().st_size:,} bytes")

    # Parse and combine
    df_all = pd.read_excel(s5_path, header=4)
    cancer_mask = df_all["Driver annotation"].str.contains("cancer", case=False, na=False)
    df_cancer = df_all[cancer_mask].copy()

    df_canonical = pd.read_excel(s3_path, header=4)
    df_canonical.columns = [
        "Gene Symbol", "Vogelstein", "Saito", "CGC",
        "Mode of action", "Cancer refs", "Alteration type", "Noncancer refs",
    ]
    role_map = dict(zip(df_canonical["Gene Symbol"], df_canonical["Mode of action"]))

    df_cancer["role"] = df_cancer["Gene Symbol"].map(role_map).fillna("unknown")
    df_cancer["is_canonical"] = df_cancer["Driver annotation"].str.contains("canonical", case=False)

    # Save as TSV
    target = output_dir / "NCG_cancerdrivers_annotation.tsv"
    df_out = df_cancer[["Gene Symbol", "Driver annotation", "role", "is_canonical", "Cancer screens (n)"]].copy()
    df_out.columns = ["ncg_symbol", "driver_annotation", "role", "is_canonical", "n_screens"]
    df_out.to_csv(target, sep="\t", index=False)
    print(f"    Built: {target} ({len(df_out)} cancer drivers)")
    return target


def download_ncg(output_dir: Path | None = None) -> Path:
    """Download NCG cancer driver gene list.

    Tries NCG website first, falls back to Springer supplementary data.

    Returns:
        Path to downloaded/built TSV file.
    """
    out = output_dir or DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    target = out / "NCG_cancerdrivers_annotation.tsv"
    if target.exists():
        # Verify it's real data, not HTML
        with open(target, "rb") as f:
            header = f.read(100)
        if not header.startswith(b"<"):
            print(f"Already exists: {target}")
            return target
        print(f"Existing file is HTML, re-downloading...")
        target.unlink()

    print("Downloading NCG cancer driver genes...")

    # Strategy 1: NCG website
    if _try_ncg_website(target):
        return target

    # Strategy 2: Springer supplementary data
    print("  NCG website unavailable, using Springer supplementary data...")
    return _build_from_springer(out)


if __name__ == "__main__":
    result = download_ncg()
    print(f"\nDone. Data at: {result}")
