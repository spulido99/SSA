"""WP8.7: Ablation study with real downloaded data.

Evaluates the impact of different configurations on driver gene prediction:
1. Baseline: mutation frequency only (2 features)
2. + Network topology (degree + betweenness)
3. + DepMap essentiality features
4. + gnomAD constraint features
5. Full model (all features)
6. Different XGBoost estimator counts (50, 100, 200)

Usage:
    python scripts/run_ablation.py
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
from ssame.data.feature_assembly import FeatureAssembler
from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.evaluation.ablation import AblationConfig, AblationResult, AblationRunner
from ssame.network.multi_view_graph import build_multi_view
from ssame.network.pioneer_oncoppi import load_pioneer


def main() -> None:
    data_dir = PROJECT_ROOT / "data" / "raw"
    t_start = time.time()

    print("=" * 60)
    print("ABLATION STUDY: SSA-ME v2 with Real Data")
    print("=" * 60)

    # ── Load mutation data ────────────────────────────────────────
    maf_path = data_dir / "tcga" / "mc3.v0.2.8.PUBLIC.maf.gz"
    if not maf_path.exists():
        print(f"SKIP: MC3 MAF not found at {maf_path}")
        return

    print("\n[1/4] Loading mutation data (first 200 samples)...")
    t0 = time.time()
    first_chunk = pd.read_csv(
        maf_path, sep="\t", comment="#",
        usecols=["Tumor_Sample_Barcode"],
        nrows=400_000,
    )
    subset_samples = set(first_chunk["Tumor_Sample_Barcode"].unique()[:200])
    alt_matrix = parse_maf_to_matrix(
        maf_path, min_samples=2, sample_barcodes=subset_samples,
    )
    print(f"  {alt_matrix.data.shape[0]} genes x {alt_matrix.data.shape[1]} samples")
    print(f"  Time: {time.time() - t0:.1f}s")

    # ── Load network ──────────────────────────────────────────────
    print("\n[2/4] Loading PIONEER network...")
    t0 = time.time()
    pioneer_path = data_dir / "networks" / "pioneer" / "pioneer_very_high_human.txt"
    if pioneer_path.exists():
        pioneer_graph = load_pioneer(pioneer_path)
        multi_view = build_multi_view({"pioneer": pioneer_graph})
        print(f"  {pioneer_graph.number_of_nodes()} nodes, {pioneer_graph.number_of_edges()} edges")
    else:
        import networkx as nx
        print("  PIONEER not found, using random graph")
        G = nx.erdos_renyi_graph(100, 0.1, seed=42)
        mapping = {i: alt_matrix.genes[i] for i in range(min(100, len(alt_matrix.genes)))}
        G = nx.relabel_nodes(G, mapping)
        multi_view = build_multi_view({"pioneer": G})
    print(f"  Time: {time.time() - t0:.1f}s")

    # ── Load ground truth ─────────────────────────────────────────
    print("\n[3/4] Building consensus ground truth...")
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
    labels_dict = consensus.as_labels()
    print(f"  {len(consensus.positives)} consensus drivers")
    print(f"  Time: {time.time() - t0:.1f}s")

    # ── Load optional features ────────────────────────────────────
    depmap_features: dict[str, float] = {}
    depmap_path = data_dir / "depmap" / "CRISPRGeneEffect.csv"
    if depmap_path.exists():
        print("  Loading DepMap essentiality...")
        df_dep = pd.read_csv(depmap_path, index_col=0)
        for col in df_dep.columns:
            gene = col.split(" (")[0] if " (" in col else col
            depmap_features[gene] = float(df_dep[col].mean())
        print(f"    {len(depmap_features)} genes with DepMap scores")

    gnomad_features: dict[str, float] = {}
    gnomad_path = data_dir / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv"
    if gnomad_path.exists():
        print("  Loading gnomAD constraint...")
        df_gn = pd.read_csv(gnomad_path, sep="\t")
        gene_col = "gene" if "gene" in df_gn.columns else df_gn.columns[0]
        loeuf_col = None
        for cand in ("oe_lof_upper", "lof.oe_ci.upper", "oe_lof_upper_bin"):
            if cand in df_gn.columns:
                loeuf_col = cand
                break
        if loeuf_col:
            for _, row in df_gn.iterrows():
                gene = str(row[gene_col]).strip()
                val = row[loeuf_col]
                if pd.notna(val):
                    gnomad_features[gene] = float(val)
        print(f"    {len(gnomad_features)} genes with gnomAD LOEUF scores")

    # ── Assemble features (different configurations) ──────────────
    print("\n[4/4] Running ablation experiments...")

    additional = {}
    if depmap_features:
        additional["depmap_essentiality"] = depmap_features
    if gnomad_features:
        additional["gnomad_loeuf"] = gnomad_features

    assembler_full = FeatureAssembler(
        alteration_matrix=alt_matrix,
        multi_view=multi_view,
        additional_features=additional,
    )
    X_full, gene_list = assembler_full.assemble(standardize=True)
    y = np.array([labels_dict.get(g, 0) for g in gene_list])

    print(f"  Feature matrix: {X_full.shape[0]} genes x {X_full.shape[1]} features")
    print(f"  Features: {assembler_full.feature_names}")
    print(f"  Positives: {int(y.sum())}/{len(y)}")

    omics_idx, topo_idx = assembler_full.get_bimodal_indices()

    # ── Feature ablations ─────────────────────────────────────────
    # Build column index map for controlled ablation
    feat_names = assembler_full.feature_names
    mut_cols = [feat_names.index("mutation_frequency"), feat_names.index("log_mutation_count")]

    ablation_runs: list[tuple[str, np.ndarray]] = []

    # 1. Mutation frequency only
    ablation_runs.append(("mut_freq_only", X_full[:, mut_cols]))

    # 2. + Network topology only
    ablation_runs.append(("+ topology", X_full[:, mut_cols + topo_idx]))

    # 3. + DepMap only (if available)
    if "depmap_essentiality" in feat_names:
        depmap_col = feat_names.index("depmap_essentiality")
        ablation_runs.append(("+ depmap", X_full[:, mut_cols + [depmap_col]]))
        ablation_runs.append(("+ depmap + topo", X_full[:, mut_cols + [depmap_col] + topo_idx]))

    # 4. + gnomAD only (if available)
    if "gnomad_loeuf" in feat_names:
        gnomad_col = feat_names.index("gnomad_loeuf")
        ablation_runs.append(("+ gnomad", X_full[:, mut_cols + [gnomad_col]]))
        ablation_runs.append(("+ gnomad + topo", X_full[:, mut_cols + [gnomad_col] + topo_idx]))

    # 5. Full model (all features)
    ablation_runs.append(("full_model", X_full))

    results: list[AblationResult] = []

    for name, X_subset in ablation_runs:
        print(f"\n  Running: {name} ({X_subset.shape[1]} features)...")
        t0 = time.time()
        runner = AblationRunner(X_subset, y)
        config = AblationConfig(
            name=name,
            encoder_level=0,
            n_estimators=100,
            n_folds=5,
        )
        result = runner.run_config(config)
        results.append(result)
        print(f"    AUPRC={result.auprc:.4f}  AUROC={result.auroc:.4f}  "
              f"P@10={result.precision_at_10:.4f}  P@50={result.precision_at_50:.4f}  "
              f"Time: {time.time() - t0:.1f}s")

    # ── Estimator count ablation ──────────────────────────────────
    for n_est in [50, 200]:
        print(f"\n  Running: full_n_est={n_est}...")
        t0 = time.time()
        runner = AblationRunner(X_full, y)
        config = AblationConfig(
            name=f"full_n_est={n_est}",
            encoder_level=0,
            n_estimators=n_est,
            n_folds=5,
        )
        result = runner.run_config(config)
        results.append(result)
        print(f"    AUPRC={result.auprc:.4f}  AUROC={result.auroc:.4f}  "
              f"P@10={result.precision_at_10:.4f}  P@50={result.precision_at_50:.4f}  "
              f"Time: {time.time() - t0:.1f}s")

    # ── Summary table ─────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("ABLATION RESULTS SUMMARY")
    print("=" * 60)
    print(AblationRunner.comparison_table(results))

    t_total = time.time() - t_start
    print(f"\nTotal time: {t_total:.1f}s")
    print("=" * 60)


if __name__ == "__main__":
    main()
