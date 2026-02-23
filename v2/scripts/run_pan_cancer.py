#!/usr/bin/env python
"""Pan-cancer per-type analysis: run SSA-ME pipeline per cancer type + Borda aggregate.

Uses TCGA barcode → cancer type mapping to split MC3 MAF by cancer type,
then runs PanCancerRunner to produce per-cancer and aggregated rankings.

Usage:
    python scripts/run_pan_cancer.py [--n-samples 100] [--n-types 10]
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.alteration_matrix import load_tcga_cancer_map, parse_maf_to_matrix
from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.network.multi_view_graph import build_multi_view
from ssame.network.pioneer_oncoppi import load_pioneer
from ssame.pipeline.pan_cancer import PanCancerRunner

RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data" / "raw"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pan-cancer analysis")
    parser.add_argument(
        "--n-samples", type=int, default=0,
        help="Max samples per cancer type (0 = all)",
    )
    parser.add_argument(
        "--n-types", type=int, default=10,
        help="Number of top cancer types to include",
    )
    parser.add_argument(
        "--n-ants", type=int, default=10,
        help="Number of ants for ACO search",
    )
    parser.add_argument(
        "--n-iter", type=int, default=50,
        help="Number of ACO iterations",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 60)
    print("PAN-CANCER PER-TYPE ANALYSIS")
    print("=" * 60)

    # ── Load barcode → cancer mapping ─────────────────────────────
    map_path = DATA_DIR / "tcga" / "tcga_barcode_cancer_map.tsv"
    if not map_path.exists():
        print(f"SKIP: Cancer map not found at {map_path}")
        print("Run: python scripts/download/download_tcga_clinical.py")
        return

    barcode_to_cancer = load_tcga_cancer_map(map_path)
    print(f"Loaded {len(barcode_to_cancer)} barcode→cancer mappings")

    # ── Load MAF barcodes and assign cancer types ─────────────────
    maf_path = DATA_DIR / "tcga" / "mc3.v0.2.8.PUBLIC.maf.gz"
    if not maf_path.exists():
        print(f"SKIP: MC3 MAF not found at {maf_path}")
        return

    print("Reading MAF sample barcodes...")
    # Read just barcodes column to determine sample→cancer mapping
    barcodes_df = pd.read_csv(
        maf_path, sep="\t", comment="#",
        usecols=["Tumor_Sample_Barcode"], nrows=5_000_000,
    )
    maf_barcodes = set(barcodes_df["Tumor_Sample_Barcode"].unique())
    print(f"  {len(maf_barcodes)} unique barcodes in MAF")

    # Map barcodes: MAF uses full barcode (TCGA-XX-XXXX-01A-01D-XXXX-XX)
    # while GDC mapping uses patient ID (TCGA-XX-XXXX).
    # Extract patient ID prefix for matching.
    cancer_barcodes: dict[str, set[str]] = {}
    unmatched = 0
    for barcode in maf_barcodes:
        # Extract patient portion: first 3 segments (TCGA-XX-XXXX)
        parts = barcode.split("-")
        if len(parts) >= 3:
            patient_id = "-".join(parts[:3])
        else:
            patient_id = barcode

        cancer_type = barcode_to_cancer.get(patient_id)
        if cancer_type is None:
            unmatched += 1
            continue

        cancer_barcodes.setdefault(cancer_type, set()).add(barcode)

    print(f"  Mapped to {len(cancer_barcodes)} cancer types ({unmatched} unmatched)")

    # Select top N cancer types by sample count
    type_counts = {ct: len(bcs) for ct, bcs in cancer_barcodes.items()}
    top_types = sorted(type_counts, key=type_counts.get, reverse=True)[:args.n_types]

    print(f"\nTop {len(top_types)} cancer types:")
    for ct in top_types:
        print(f"  {ct}: {type_counts[ct]} samples")

    # ── Load ground truth ─────────────────────────────────────────
    gt_dir = DATA_DIR / "ground_truth"
    intogen_path = gt_dir / "intogen" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"
    ncg_path = gt_dir / "NCG_cancerdrivers_annotation.tsv"
    cm_path = gt_dir / "cancermine_collated.tsv"

    consensus = build_consensus_ground_truth(
        intogen_path=intogen_path if intogen_path.exists() else None,
        ncg_path=ncg_path if ncg_path.exists() else None,
        cancermine_path=cm_path if cm_path.exists() else None,
        min_sources=2,
    )
    labels = consensus.as_labels()
    print(f"Ground truth: {len(labels)} consensus drivers")

    # ── Parse MAF per cancer type ─────────────────────────────────
    from ssame.types import AlterationMatrix

    alt_matrices: dict[str, AlterationMatrix] = {}
    for cancer_type in top_types:
        barcodes = cancer_barcodes[cancer_type]
        if args.n_samples > 0 and len(barcodes) > args.n_samples:
            barcodes = set(sorted(barcodes)[:args.n_samples])

        print(f"\nParsing MAF for {cancer_type} ({len(barcodes)} samples)...")
        matrix = parse_maf_to_matrix(
            maf_path,
            min_samples=2,
            sample_barcodes=barcodes,
        )
        if matrix.n_genes > 0 and matrix.n_samples > 0:
            alt_matrices[cancer_type] = matrix
            print(f"  {matrix.n_genes} genes x {matrix.n_samples} samples")
        else:
            print(f"  SKIP: empty matrix")

    if not alt_matrices:
        print("ERROR: No valid cancer-type matrices produced.")
        return

    # ── Run pan-cancer ────────────────────────────────────────────
    print(f"\nRunning PanCancerRunner on {len(alt_matrices)} cancer types...")
    runner = PanCancerRunner(
        alteration_matrices=alt_matrices,
        labels=labels,
        encoder_level=0,
        encoder_n_estimators=100,
        search_n_iterations=args.n_iter,
        seed=42,
    )
    result = runner.run()

    # ── Save results ──────────────────────────────────────────────
    # Aggregated ranking
    rows = []
    for rank, (gene, score) in enumerate(result.aggregated_ranking, 1):
        rows.append({
            "rank": rank,
            "gene": gene,
            "borda_score": f"{score:.6f}",
            "is_known_driver": gene in labels,
            "recurrence": result.recurrence_count.get(gene, 0),
        })
    df_agg = pd.DataFrame(rows)
    df_agg.to_csv(RESULTS_DIR / "pan_cancer_aggregated.tsv", sep="\t", index=False)
    print(f"\nSaved: pan_cancer_aggregated.tsv ({len(rows)} genes)")

    # Per-cancer type top genes
    per_type_rows = []
    for cancer_type, pr in result.per_cancer.items():
        for rank, (gene, score) in enumerate(pr.ranked_genes[:100], 1):
            per_type_rows.append({
                "cancer_type": cancer_type,
                "rank": rank,
                "gene": gene,
                "score": f"{score:.6f}",
                "is_known_driver": gene in labels,
            })
    df_per = pd.DataFrame(per_type_rows)
    df_per.to_csv(RESULTS_DIR / "pan_cancer_per_type.tsv", sep="\t", index=False)
    print(f"Saved: pan_cancer_per_type.tsv")

    # Recurrence
    rec_rows = []
    for gene, count in sorted(result.recurrence_count.items(), key=lambda x: -x[1]):
        if count > 0:
            rec_rows.append({
                "gene": gene,
                "n_cancer_types": count,
                "is_known_driver": gene in labels,
            })
    df_rec = pd.DataFrame(rec_rows)
    df_rec.to_csv(RESULTS_DIR / "pan_cancer_recurrence.tsv", sep="\t", index=False)
    print(f"Saved: pan_cancer_recurrence.tsv ({len(rec_rows)} recurrent genes)")

    # Summary stats
    top10_agg = result.aggregated_ranking[:10]
    n_drivers_top10 = sum(1 for g, _ in top10_agg if g in labels)
    top50_agg = result.aggregated_ranking[:50]
    n_drivers_top50 = sum(1 for g, _ in top50_agg if g in labels)

    print(f"\n--- Summary ---")
    print(f"Cancer types analyzed: {len(result.per_cancer)}")
    print(f"Aggregated top-10: {n_drivers_top10}/10 known drivers")
    print(f"Aggregated top-50: {n_drivers_top50}/50 known drivers")
    print(f"Genes recurring in >=3 types: {sum(1 for c in result.recurrence_count.values() if c >= 3)}")

    t_total = time.time() - t_start
    print(f"\nPan-cancer analysis completed in {t_total:.1f}s")


if __name__ == "__main__":
    main()
