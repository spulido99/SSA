#!/usr/bin/env python
"""Bias-stratified evaluation: assess whether SSA-ME predictions are biased toward well-studied genes.

Uses gnomAD LOEUF as study-bias proxy (lower LOEUF = more constrained = more studied)
and mutation frequency as secondary proxy, then runs BiasStratifiedEvaluator
with 4 quartiles to measure precision across strata.

Usage:
    python scripts/run_bias_evaluation.py
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.evaluation.bias_stratified import BiasStratifiedEvaluator

RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data" / "raw"


def load_ranked_genes(path: Path) -> dict[str, float]:
    """Load ranked_genes.tsv -> dict[gene, score]."""
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["gene"], df["score"].astype(float)))


def load_gnomad_loeuf(path: Path) -> dict[str, float]:
    """Load gnomAD LOEUF scores as study-bias proxy.

    Lower LOEUF = more constrained = typically more studied.
    We invert to: higher value = more well-studied.
    """
    df = pd.read_csv(path, sep="\t")
    loeuf: dict[str, float] = {}
    for _, row in df.iterrows():
        gene = str(row["gene"]).strip()
        val = row.get("lof.oe_ci.upper")
        if pd.notna(val):
            # Invert: well-studied genes have LOW LOEUF, so use 1/LOEUF
            fval = float(val)
            if fval > 0:
                loeuf[gene] = 1.0 / fval
    return loeuf


def load_mutation_frequency(path: Path) -> dict[str, float]:
    """Load mutation frequency from ranked_genes.tsv as alternative bias proxy.

    Genes with higher mutation frequency tend to be more well-studied.
    This is computed from the alteration matrix in the pipeline.
    """
    df = pd.read_csv(path, sep="\t")
    # Use rank as proxy (lower rank = higher frequency = more studied)
    # Invert: rank 1 gets highest value
    n = len(df)
    return dict(zip(df["gene"], [n - r for r in range(n)]))


def run_single_evaluation(
    gene_scores: dict[str, float],
    study_bias: dict[str, float],
    labels: dict[str, int],
    proxy_name: str,
    n_strata: int = 4,
) -> dict:
    """Run a single bias-stratified evaluation and return results."""
    evaluator = BiasStratifiedEvaluator(
        gene_scores=gene_scores,
        study_bias_proxy=study_bias,
        labels=labels,
        n_strata=n_strata,
    )
    result = evaluator.evaluate()

    print(f"\n--- {proxy_name} ---")
    print(f"Score-bias correlation: {result.summary['score_bias_correlation']:.4f} "
          f"(p={result.summary['score_bias_pvalue']:.2e})")
    print(f"Genes evaluated: {result.summary['n_genes']}")

    print(f"\n{'Stratum':>8} {'N genes':>8} {'N pos':>6} {'Mean bias':>10} "
          f"{'Mean score':>11} {'P@k':>6}")
    print("-" * 55)

    strata_rows = []
    for s in result.strata_metrics:
        print(f"{s['stratum']:>8} {s['n_genes']:>8} {s['n_positives']:>6} "
              f"{s['mean_study_bias']:>10.1f} {s['mean_score']:>11.4f} "
              f"{s['precision_at_k']:>6.3f}")
        strata_rows.append({
            "proxy": proxy_name,
            "stratum": s["stratum"],
            "n_genes": s["n_genes"],
            "n_positives": s["n_positives"],
            "mean_study_bias": f"{s['mean_study_bias']:.2f}",
            "mean_score": f"{s['mean_score']:.6f}",
            "precision_at_k": f"{s['precision_at_k']:.4f}",
        })

    return {
        "proxy_name": proxy_name,
        "result": result,
        "strata_rows": strata_rows,
    }


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

    all_strata_rows = []

    # ── Proxy 1: gnomAD LOEUF (inverse = constraint) ─────────────
    gnomad_path = DATA_DIR / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv"
    if gnomad_path.exists():
        print("\nLoading gnomAD LOEUF as study-bias proxy...")
        loeuf_proxy = load_gnomad_loeuf(gnomad_path)
        print(f"  {len(loeuf_proxy)} genes with LOEUF scores")

        overlap = set(gene_scores) & set(loeuf_proxy)
        print(f"  Overlap with gene scores: {len(overlap)} genes")

        if len(overlap) > 10:
            res = run_single_evaluation(
                gene_scores, loeuf_proxy, labels, "gnomAD_LOEUF_inverse"
            )
            all_strata_rows.extend(res["strata_rows"])
    else:
        print(f"\nSKIP gnomAD: {gnomad_path} not found")

    # ── Proxy 2: Mutation frequency rank ──────────────────────────
    print("\nUsing mutation frequency rank as study-bias proxy...")
    freq_proxy = load_mutation_frequency(ranked_path)
    print(f"  {len(freq_proxy)} genes with frequency proxy")

    res2 = run_single_evaluation(
        gene_scores, freq_proxy, labels, "mutation_frequency_rank"
    )
    all_strata_rows.extend(res2["strata_rows"])

    # ── Save all results ──────────────────────────────────────────
    df_strata = pd.DataFrame(all_strata_rows)
    df_strata.to_csv(RESULTS_DIR / "bias_stratified.tsv", sep="\t", index=False)
    print(f"\nSaved: bias_stratified.tsv")

    # Save summary
    summary_items = {}
    for res_data in [res.get("result") for res in [res2] if "result" in res]:
        pass
    # Use the last available result for summary
    last_result = res2["result"]
    summary_items = {
        "score_bias_correlation": f"{last_result.summary['score_bias_correlation']:.6f}",
        "score_bias_pvalue": f"{last_result.summary['score_bias_pvalue']:.2e}",
        "n_genes": str(int(last_result.summary["n_genes"])),
        "n_strata": str(int(last_result.summary["n_strata"])),
        "bias_proxies_used": "gnomAD_LOEUF_inverse,mutation_frequency_rank",
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
