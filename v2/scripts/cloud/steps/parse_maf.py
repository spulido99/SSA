#!/usr/bin/env python
"""SageMaker processing step: Parse MAF into per-cancer-type alteration matrices.

Reads MAF from /opt/ml/processing/input/
Writes per-cancer alteration matrices to /opt/ml/processing/output/{cancer_type}/
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import joblib
import pandas as pd

# Ensure ssame package is importable
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.alteration_matrix import parse_maf_to_matrix

INPUT_DIR = Path(os.environ.get("SM_INPUT_DIR", "/opt/ml/processing/input"))
OUTPUT_DIR = Path(os.environ.get("SM_OUTPUT_DIR", "/opt/ml/processing/output"))


def main() -> None:
    t_start = time.time()
    print("=" * 60)
    print("STEP 1: ParseMAF")
    print("=" * 60)

    # Find MAF file
    maf_candidates = list(INPUT_DIR.rglob("*.maf")) + list(INPUT_DIR.rglob("*.maf.gz"))
    if not maf_candidates:
        print(f"ERROR: No MAF file found in {INPUT_DIR}")
        sys.exit(1)
    maf_path = maf_candidates[0]
    print(f"Using MAF: {maf_path}")

    # Find barcode master
    master_candidates = list(INPUT_DIR.rglob("barcode_master.tsv"))
    if not master_candidates:
        print(f"ERROR: No barcode_master.tsv found in {INPUT_DIR}")
        sys.exit(1)
    master_df = pd.read_csv(master_candidates[0], sep="\t")
    print(f"Loaded barcode master: {len(master_df)} barcodes")

    # Group by cancer type
    cancer_barcodes = master_df.groupby("cancer_type")["barcode"].apply(set).to_dict()

    n_types = int(os.environ.get("N_CANCER_TYPES", "10"))
    type_counts = {ct: len(bcs) for ct, bcs in cancer_barcodes.items()}
    top_types = sorted(type_counts, key=type_counts.get, reverse=True)[:n_types]

    print(f"\nParsing {len(top_types)} cancer types...")
    for cancer_type in top_types:
        barcodes = cancer_barcodes[cancer_type]
        t_parse = time.time()
        print(f"  {cancer_type} ({len(barcodes)} samples)...", end=" ", flush=True)

        matrix = parse_maf_to_matrix(maf_path, min_samples=2, sample_barcodes=barcodes)

        if matrix.n_genes > 0 and matrix.n_samples > 0:
            out_dir = OUTPUT_DIR / cancer_type
            out_dir.mkdir(parents=True, exist_ok=True)
            joblib.dump(matrix, out_dir / "alteration_matrix.pkl")
            dt = time.time() - t_parse
            print(f"{matrix.n_genes} genes x {matrix.n_samples} samples ({dt:.1f}s)")
        else:
            print("SKIP: empty matrix")

    dt_total = time.time() - t_start
    print(f"\nParseMAF completed in {dt_total:.1f}s")


if __name__ == "__main__":
    main()
