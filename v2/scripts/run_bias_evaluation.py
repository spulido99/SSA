#!/usr/bin/env python
"""Bias-stratified evaluation: assess whether SSA-ME predictions are biased toward well-studied genes.

Uses PIONEER network degree as study-bias proxy and runs BiasStratifiedEvaluator
with 4 quartiles to measure precision across strata.

Usage:
    python scripts/run_bias_evaluation.py
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import networkx as nx
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.evaluation.bias_stratified import BiasStratifiedEvaluator
from ssame.network.pioneer_oncoppi import load_pioneer

RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data" / "raw"


def load_ranked_genes(path: Path) -> dict[str, float]:
    """Load ranked_genes.tsv → dict[gene, score]."""
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["gene"], df["score"].astype(float)))


def compute_network_degree(graph: nx.Graph) -> dict[str, int]:
    """Compute node degree in network as study-bias proxy."""
    return dict(graph.degree())


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 60)
    print("BIAS-STRATIFIED EVALUATION")
    print("=" * 60)

    # ── Load gene scores ──────────────────────────────────────────
    ranked_path = RESULTS_DIR / "ranked_genes.tsv"
    if not ranked_path.exists():
        print(f"SKIP: ranked_genes.tsv not found at {ranked_path}")
        print("Run scripts/save_results.py first.")
        return

    gene_scores = load_ranked_genes(ranked_path)
    print(f"Loaded {len(gene_scores)} gene scores")

    # ── Load PIONEER network for study-bias proxy ─────────────────
    pioneer_path = DATA_DIR / "networks" / "pioneer" / "pioneer_very_high_human.txt"
    if not pioneer_path.exists():
        print(f"SKIP: PIONEER network not found at {pioneer_path}")
        return

    print("Loading PIONEER network...")
    pioneer_graph = load_pioneer(pioneer_path)
    network_degree = compute_network_degree(pioneer_graph)
    print(f"  {len(network_degree)} genes with degree info")

    # ── Load ground truth labels ──────────────────────────────────
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

    # ── Run bias-stratified evaluation ────────────────────────────
    print("\nRunning bias-stratified evaluation (4 quartiles)...")
    evaluator = BiasStratifiedEvaluator(
        gene_scores=gene_scores,
        study_bias_proxy=network_degree,
        labels=labels,
        n_strata=4,
    )
    result = evaluator.evaluate()

    # ── Display results ───────────────────────────────────────────
    print(f"\nScore-bias correlation: {result.summary['score_bias_correlation']:.4f} "
          f"(p={result.summary['score_bias_pvalue']:.2e})")
    print(f"Total genes evaluated: {result.summary['n_genes']}")

    print("\nPer-stratum metrics:")
    print(f"{'Stratum':>8} {'N genes':>8} {'N pos':>6} {'Mean bias':>10} "
          f"{'Mean score':>11} {'P@k':>6}")
    print("-" * 55)

    strata_rows = []
    for s in result.strata_metrics:
        print(f"{s['stratum']:>8} {s['n_genes']:>8} {s['n_positives']:>6} "
              f"{s['mean_study_bias']:>10.1f} {s['mean_score']:>11.4f} "
              f"{s['precision_at_k']:>6.3f}")
        strata_rows.append({
            "stratum": s["stratum"],
            "n_genes": s["n_genes"],
            "n_positives": s["n_positives"],
            "mean_study_bias": f"{s['mean_study_bias']:.2f}",
            "mean_score": f"{s['mean_score']:.6f}",
            "precision_at_k": f"{s['precision_at_k']:.4f}",
        })

    # ── Save results ──────────────────────────────────────────────
    df_strata = pd.DataFrame(strata_rows)
    df_strata.to_csv(RESULTS_DIR / "bias_stratified.tsv", sep="\t", index=False)
    print(f"\nSaved: bias_stratified.tsv")

    # Save summary
    summary_items = {
        "score_bias_correlation": f"{result.summary['score_bias_correlation']:.6f}",
        "score_bias_pvalue": f"{result.summary['score_bias_pvalue']:.2e}",
        "n_genes": str(int(result.summary["n_genes"])),
        "n_strata": str(int(result.summary["n_strata"])),
        "bias_proxy": "PIONEER_network_degree",
    }
    with open(RESULTS_DIR / "bias_summary.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        for k, v in summary_items.items():
            writer.writerow([k, v])
    print(f"Saved: bias_summary.tsv")

    t_total = time.time() - t_start
    print(f"\nBias evaluation completed in {t_total:.1f}s")


if __name__ == "__main__":
    main()
