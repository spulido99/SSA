#!/usr/bin/env python
"""SageMaker processing step: Aggregate per-cancer results via Borda ranking.

Reads per-cancer PipelineResult files from /opt/ml/processing/input/
Writes final TSVs to /opt/ml/processing/output/
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

from ssame.pipeline.pan_cancer import PanCancerRunner
from ssame.pipeline.pipeline import PipelineResult

INPUT_DIR = Path(os.environ.get("SM_INPUT_DIR", "/opt/ml/processing/input"))
OUTPUT_DIR = Path(os.environ.get("SM_OUTPUT_DIR", "/opt/ml/processing/output"))


def main() -> None:
    t_start = time.time()
    print("=" * 60)
    print("STEP 3: Aggregate")
    print("=" * 60)

    # Load all per-cancer results
    per_cancer: dict[str, PipelineResult] = {}
    for result_file in sorted(INPUT_DIR.rglob("pipeline_result.pkl")):
        cancer_type = result_file.parent.name
        per_cancer[cancer_type] = joblib.load(result_file)
        n_genes = len(per_cancer[cancer_type].ranked_genes)
        print(f"  Loaded {cancer_type}: {n_genes} genes")

    if not per_cancer:
        print("ERROR: No per-cancer results found")
        sys.exit(1)

    print(f"\nAggregating {len(per_cancer)} cancer types...")

    # Use PanCancerRunner's aggregation methods (create a dummy runner)
    runner = PanCancerRunner.__new__(PanCancerRunner)
    runner._top_k_fraction = 0.1
    aggregated = runner._borda_aggregate(per_cancer)
    recurrence = runner._compute_recurrence(per_cancer)

    # Save results
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Aggregated ranking
    rows = []
    for rank, (gene, score) in enumerate(aggregated, 1):
        rows.append({
            "rank": rank,
            "gene": gene,
            "borda_score": f"{score:.6f}",
            "recurrence": recurrence.get(gene, 0),
        })
    df_agg = pd.DataFrame(rows)
    df_agg.to_csv(OUTPUT_DIR / "pan_cancer_aggregated.tsv", sep="\t", index=False)
    print(f"Saved: pan_cancer_aggregated.tsv ({len(rows)} genes)")

    # Per-cancer type top genes
    per_type_rows = []
    for cancer_type, pr in per_cancer.items():
        for rank, (gene, score) in enumerate(pr.ranked_genes[:100], 1):
            per_type_rows.append({
                "cancer_type": cancer_type,
                "rank": rank,
                "gene": gene,
                "score": f"{score:.6f}",
            })
    df_per = pd.DataFrame(per_type_rows)
    df_per.to_csv(OUTPUT_DIR / "pan_cancer_per_type.tsv", sep="\t", index=False)
    print(f"Saved: pan_cancer_per_type.tsv")

    # Recurrence
    rec_rows = []
    for gene, count in sorted(recurrence.items(), key=lambda x: -x[1]):
        if count > 0:
            rec_rows.append({"gene": gene, "n_cancer_types": count})
    df_rec = pd.DataFrame(rec_rows)
    df_rec.to_csv(OUTPUT_DIR / "pan_cancer_recurrence.tsv", sep="\t", index=False)
    print(f"Saved: pan_cancer_recurrence.tsv ({len(rec_rows)} recurrent genes)")

    dt = time.time() - t_start
    print(f"\nAggregation completed in {dt:.1f}s")


if __name__ == "__main__":
    main()
