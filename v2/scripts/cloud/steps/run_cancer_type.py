#!/usr/bin/env python
"""SageMaker processing step: Run SSA-ME pipeline for one cancer type.

Reads alteration matrix from /opt/ml/processing/input/
Writes PipelineResult to /opt/ml/processing/output/
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import joblib

# Ensure ssame package is importable
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.pipeline.pipeline import Pipeline
from ssame.pipeline.pan_cancer import PanCancerRunner

INPUT_DIR = Path(os.environ.get("SM_INPUT_DIR", "/opt/ml/processing/input"))
OUTPUT_DIR = Path(os.environ.get("SM_OUTPUT_DIR", "/opt/ml/processing/output"))


def main() -> None:
    cancer_type = os.environ.get("CANCER_TYPE", "UNKNOWN")
    encoder_level = int(os.environ.get("ENCODER_LEVEL", "0"))
    n_iter = int(os.environ.get("N_ITER", "200"))
    n_ants = int(os.environ.get("N_ANTS", "20"))

    t_start = time.time()
    print(f"Running pipeline for {cancer_type} (L{encoder_level})")

    # Load alteration matrix
    alt_path = INPUT_DIR / "alteration_matrix.pkl"
    if not alt_path.exists():
        print(f"ERROR: Alteration matrix not found at {alt_path}")
        sys.exit(1)
    alt_matrix = joblib.load(alt_path)
    print(f"  Loaded: {alt_matrix.n_genes} genes x {alt_matrix.n_samples} samples")

    # Build ground truth labels
    gt_dir = INPUT_DIR / "ground_truth"
    if gt_dir.exists():
        intogen_path = next(gt_dir.rglob("Compendium_Cancer_Genes.tsv"), None)
        ncg_path = next(gt_dir.rglob("NCG_cancerdrivers_annotation.tsv"), None)
        cm_path = next(gt_dir.rglob("cancermine_collated.tsv"), None)
        consensus = build_consensus_ground_truth(
            intogen_path=intogen_path,
            ncg_path=ncg_path,
            cancermine_path=cm_path,
            min_sources=2,
        )
        labels = consensus.as_labels()
    else:
        print("  WARNING: No ground truth directory found, using empty labels")
        labels = {}

    # Run pipeline with checkpoint support
    ckpt_dir = OUTPUT_DIR / ".checkpoints"
    runner = PanCancerRunner(
        alteration_matrices={cancer_type: alt_matrix},
        labels=labels,
        encoder_level=encoder_level,
        search_n_iterations=n_iter,
        seed=42,
        checkpoint_dir=ckpt_dir,
    )
    result = runner.run()

    # Save result
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cancer_result = result.per_cancer[cancer_type]
    joblib.dump(cancer_result, OUTPUT_DIR / "pipeline_result.pkl")

    dt = time.time() - t_start
    n_genes = len(cancer_result.ranked_genes)
    top3 = [g for g, _ in cancer_result.ranked_genes[:3]]
    print(f"  Done: {n_genes} genes ranked, top 3: {top3} ({dt:.1f}s)")


if __name__ == "__main__":
    main()
