#!/usr/bin/env python
"""Pan-cancer per-type analysis: run SSA-ME pipeline per cancer type + Borda aggregate.

Uses barcode master cache for cancer-type splitting and uncompressed MAF
with PyArrow for fast per-type parsing.

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

from ssame.data.alteration_matrix import parse_maf_to_matrix
from ssame.data.open_ground_truth import build_consensus_ground_truth
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
    parser.add_argument("--n-workers", type=int, default=0, help="Thread limit (0=auto)")
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from checkpoints (skip completed cancer types).",
    )
    parser.add_argument(
        "--fresh", action="store_true",
        help="Clear all checkpoints before running.",
    )
    parser.add_argument(
        "--mode", choices=["lite", "standard", "full"], default=None,
        help="Apply parameter preset (overrides --n-ants/--n-iter).",
    )
    return parser.parse_args()


def load_barcode_master(path: Path) -> pd.DataFrame:
    """Load barcode master cache (barcode, patient_id, cancer_type)."""
    if not path.exists():
        raise FileNotFoundError(
            f"Barcode master not found: {path}\n"
            "Run: python -c \"...\" to generate it (see run_pan_cancer.py header)"
        )
    return pd.read_csv(path, sep="\t")


def main() -> None:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    # Apply preset if given
    from ssame.config import PRESETS
    if args.mode:
        preset = PRESETS[args.mode]
        args.n_ants = preset.get("n_ants", args.n_ants)
        args.n_iter = preset.get("n_iter", args.n_iter)

    print("=" * 60)
    print("PAN-CANCER PER-TYPE ANALYSIS")
    print("=" * 60)

    # ── Load barcode master cache ─────────────────────────────────
    master_path = DATA_DIR / "tcga" / "barcode_master.tsv"
    if not master_path.exists():
        print(f"SKIP: Barcode master not found at {master_path}")
        print("Generate it first (see docs).")
        return

    master_df = load_barcode_master(master_path)
    print(f"Loaded barcode master: {len(master_df)} barcodes")

    # Group barcodes by cancer type (vectorized)
    cancer_barcodes = master_df.groupby("cancer_type")["barcode"].apply(set).to_dict()

    # Select top N cancer types by sample count
    type_counts = {ct: len(bcs) for ct, bcs in cancer_barcodes.items()}
    top_types = sorted(type_counts, key=type_counts.get, reverse=True)[:args.n_types]

    print(f"\nTop {len(top_types)} cancer types:")
    for ct in top_types:
        print(f"  {ct}: {type_counts[ct]} samples")

    # ── Resolve MAF path (prefer uncompressed for PyArrow) ────────
    maf_uncompressed = DATA_DIR / "tcga" / "mc3.v0.2.8.PUBLIC.maf"
    maf_compressed = DATA_DIR / "tcga" / "mc3.v0.2.8.PUBLIC.maf.gz"

    if maf_uncompressed.exists():
        maf_path = maf_uncompressed
        print(f"\nUsing uncompressed MAF (PyArrow fast path)")
    elif maf_compressed.exists():
        maf_path = maf_compressed
        print(f"\nUsing compressed MAF (pandas, slower)")
    else:
        print("SKIP: MC3 MAF not found")
        return

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

        t_parse = time.time()
        print(f"\nParsing MAF for {cancer_type} ({len(barcodes)} samples)...", end=" ", flush=True)
        matrix = parse_maf_to_matrix(
            maf_path,
            min_samples=2,
            sample_barcodes=barcodes,
        )
        dt = time.time() - t_parse
        if matrix.n_genes > 0 and matrix.n_samples > 0:
            alt_matrices[cancer_type] = matrix
            print(f"{matrix.n_genes} genes x {matrix.n_samples} samples ({dt:.1f}s)")
        else:
            print(f"SKIP: empty matrix ({dt:.1f}s)")

    if not alt_matrices:
        print("ERROR: No valid cancer-type matrices produced.")
        return

    # ── Checkpoint setup ──────────────────────────────────────────
    checkpoint_dir = None
    if args.resume:
        checkpoint_dir = RESULTS_DIR / ".checkpoints_pan_cancer"
    if args.fresh and checkpoint_dir and checkpoint_dir.exists():
        import shutil
        shutil.rmtree(checkpoint_dir)
        print("Cleared all checkpoints")

    # ── Run pan-cancer ────────────────────────────────────────────
    print(f"\nRunning PanCancerRunner on {len(alt_matrices)} cancer types...")
    runner = PanCancerRunner(
        alteration_matrices=alt_matrices,
        labels=labels,
        encoder_level=0,
        encoder_n_estimators=100,
        search_n_iterations=args.n_iter,
        seed=42,
        checkpoint_dir=checkpoint_dir,
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
