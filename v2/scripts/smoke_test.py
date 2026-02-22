"""Smoke test: run the full pipeline with real downloaded data.

Uses a small subset (100 samples) for speed. Validates:
1. MAF parser works with compressed MC3 data
2. Network loaders work with PIONEER data
3. Ground truth consensus builds correctly
4. Pipeline L0 (XGBoost) runs end-to-end
5. Output contains ranked genes and modules

Usage:
    python scripts/smoke_test.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from ssame.data.alteration_matrix import parse_maf_to_matrix
from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.network.multi_view_graph import build_multi_view
from ssame.network.pioneer_oncoppi import load_pioneer
from ssame.pipeline.pipeline import Pipeline


def main() -> None:
    data_dir = PROJECT_ROOT / "data" / "raw"
    t_start = time.time()

    # ── Step 1: Select a subset of samples ──────────────────────────
    print("=" * 60)
    print("SMOKE TEST: SSA-ME v2 Pipeline with Real Data")
    print("=" * 60)

    maf_path = data_dir / "tcga" / "mc3.v0.2.8.PUBLIC.maf.gz"
    if not maf_path.exists():
        print(f"SKIP: MC3 MAF not found at {maf_path}")
        return

    print("\n[1/5] Loading mutation data (first 100 samples)...")
    t0 = time.time()

    # Read first chunk to get 100 unique sample barcodes
    first_chunk = pd.read_csv(
        maf_path, sep="\t", comment="#",
        usecols=["Tumor_Sample_Barcode"],
        nrows=200_000,
    )
    subset_samples = set(first_chunk["Tumor_Sample_Barcode"].unique()[:100])
    print(f"  Selected {len(subset_samples)} samples")

    # Parse MAF for subset only
    alt_matrix = parse_maf_to_matrix(
        maf_path,
        min_samples=2,
        sample_barcodes=subset_samples,
    )
    t1 = time.time()
    print(f"  AlterationMatrix: {alt_matrix.data.shape[0]} genes x {alt_matrix.data.shape[1]} samples")
    print(f"  Density: {alt_matrix.data.sum() / alt_matrix.data.size:.4f}")
    print(f"  Time: {t1 - t0:.1f}s")

    # ── Step 2: Build network ────────────────────────────────────────
    print("\n[2/5] Building multi-view network...")
    t0 = time.time()

    pioneer_path = data_dir / "networks" / "pioneer" / "pioneer_very_high_human.txt"
    if pioneer_path.exists():
        pioneer_graph = load_pioneer(pioneer_path)
        print(f"  PIONEER: {pioneer_graph.number_of_nodes()} nodes, {pioneer_graph.number_of_edges()} edges")
    else:
        import networkx as nx
        print("  PIONEER not found, using random graph")
        pioneer_graph = nx.erdos_renyi_graph(100, 0.1, seed=42)
        mapping = {i: alt_matrix.genes[i] for i in range(min(100, len(alt_matrix.genes)))}
        pioneer_graph = nx.relabel_nodes(pioneer_graph, mapping)

    multi_view = build_multi_view({"pioneer": pioneer_graph})
    t1 = time.time()
    print(f"  MultiViewGraph: {len(multi_view.node_list)} nodes, views={multi_view.view_names}")
    print(f"  Time: {t1 - t0:.1f}s")

    # ── Step 3: Build ground truth labels ────────────────────────────
    print("\n[3/5] Building consensus ground truth...")
    t0 = time.time()

    gt_dir = data_dir / "ground_truth"
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
    t1 = time.time()
    print(f"  Consensus: {len(consensus.positives)} drivers (>=2 sources)")
    print(f"  High confidence: {len(consensus.high_confidence)} (all 3 sources)")

    # How many labels overlap with alteration matrix genes?
    labeled_in_matrix = sum(1 for g in alt_matrix.genes if g in labels)
    print(f"  Labeled genes in matrix: {labeled_in_matrix}/{alt_matrix.data.shape[0]}")
    print(f"  Time: {t1 - t0:.1f}s")

    # ── Step 4: Run Pipeline (L0 = XGBoost) ──────────────────────────
    print("\n[4/5] Running Pipeline (encoder_level=0, XGBoost)...")
    t0 = time.time()

    pipeline = Pipeline(
        alteration_matrix=alt_matrix,
        multi_view=multi_view,
        labels=labels,
        encoder_level=0,
        encoder_n_estimators=50,
        search_n_ants=10,
        search_n_iterations=50,
        search_subnetwork_size=3,
        seed=42,
    )
    result = pipeline.run()
    t1 = time.time()
    print(f"  Ranked genes: {len(result.ranked_genes)}")
    print(f"  Modules found: {len(result.modules)}")
    print(f"  Module scores: {len(result.module_scores)}")
    print(f"  Time: {t1 - t0:.1f}s")

    # ── Step 5: Validate output ──────────────────────────────────────
    print("\n[5/5] Validating output...")

    # Check known drivers are ranked high
    known_in_top100 = sum(
        1 for gene, _score in result.ranked_genes[:100]
        if gene in labels
    )
    total_known = sum(1 for gene, _ in result.ranked_genes if gene in labels)
    print(f"  Known drivers in top 100: {known_in_top100}/{min(100, len(result.ranked_genes))}")
    print(f"  Total known drivers ranked: {total_known}")

    # Top 20 genes
    print("\n  Top 20 ranked genes:")
    for i, (gene, score) in enumerate(result.ranked_genes[:20]):
        marker = "*" if gene in labels else " "
        print(f"    {i + 1:3d}. {marker} {gene:12s} score={score:.4f}")

    print(f"\n  (* = known driver in consensus ground truth)")

    # Sanity checks
    assert len(result.ranked_genes) > 0, "No ranked genes produced"
    assert len(result.modules) > 0, "No modules found"
    assert all(isinstance(s, float) for s in result.module_scores), "Invalid module scores"

    # Check scores are reasonable
    scores = [s for _, s in result.ranked_genes]
    assert not np.isnan(scores).any(), "NaN in gene scores"
    assert max(scores) > min(scores), "All scores identical"

    t_total = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"SMOKE TEST PASSED in {t_total:.1f}s")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
