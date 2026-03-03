#!/usr/bin/env python
"""Per-cancer-type ablation: show mechanistic features matter most in low-mutation cancers.

For each cancer type, compares feature configurations:
1. mut_only: mutation frequency alone
2. + gnomad: adds gnomAD LOEUF constraint
3. + depmap: adds DepMap CRISPR essentiality
4. + topology: adds co-alteration network degree/betweenness
5. full: all features combined

Hypothesis: in low-mutation cancers (THCA, PRAD, LGG), mechanistic features
(gnomAD, DepMap) provide the largest relative lift over mutation-only baseline,
demonstrating the method's power beyond simple recurrence counting.

Usage:
    python scripts/run_per_type_ablation.py [--n-types 10]
"""

from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scipy.sparse import csr_matrix

from ssame.data.alteration_matrix import parse_maf_to_matrix
from ssame.data.feature_assembly import FeatureAssembler
from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.evaluation.ablation import AblationConfig, AblationResult, AblationRunner
from ssame.network.multi_view_graph import MultiViewGraph

RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data" / "raw"


def build_coalteration_multi_view(
    alt_matrix_data: np.ndarray, genes: list[str], n_samples: int,
) -> MultiViewGraph:
    """Build co-alteration multi-view graph for feature assembly."""
    threshold = max(2, n_samples // 50)
    sparse_data = csr_matrix(alt_matrix_data.astype(np.float32))

    if len(genes) <= 8000:
        co = sparse_data @ sparse_data.T
        co.setdiag(0)
        co.eliminate_zeros()
        adj = (co >= threshold).astype(np.float32)
        adj.eliminate_zeros()
        del co
    else:
        from scipy.sparse import vstack
        chunk_size = 2000
        blocks = []
        for start in range(0, len(genes), chunk_size):
            end = min(start + chunk_size, len(genes))
            block = sparse_data[start:end] @ sparse_data.T
            block = (block >= threshold).astype(np.float32)
            block.eliminate_zeros()
            blocks.append(block.tocsr())
            del block
        adj = vstack(blocks, format="csr")
        adj.setdiag(0)
        adj.eliminate_zeros()

    node_to_index = {g: i for i, g in enumerate(genes)}
    return MultiViewGraph(
        node_list=genes,
        node_to_index=node_to_index,
        view_names=["coalteration"],
        _adjacencies={"coalteration": adj},
    )


def run_ablations_for_type(
    cancer_type: str,
    X_full: np.ndarray,
    y: np.ndarray,
    feature_names: list[str],
    omics_idx: list[int],
    topo_idx: list[int],
) -> list[AblationResult]:
    """Run feature ablation experiments for a single cancer type."""
    mut_cols = [feature_names.index("mutation_frequency"), feature_names.index("log_mutation_count")]

    configs: list[tuple[str, np.ndarray]] = []

    # 1. Mutation only
    configs.append(("mut_only", X_full[:, mut_cols]))

    # 2. + gnomAD
    if "gnomad_loeuf" in feature_names:
        gnomad_col = feature_names.index("gnomad_loeuf")
        configs.append(("+ gnomad", X_full[:, mut_cols + [gnomad_col]]))

    # 3. + DepMap
    if "depmap_essentiality" in feature_names:
        depmap_col = feature_names.index("depmap_essentiality")
        configs.append(("+ depmap", X_full[:, mut_cols + [depmap_col]]))

    # 4. + topology
    if topo_idx:
        configs.append(("+ topology", X_full[:, mut_cols + topo_idx]))

    # 5. Full
    configs.append(("full", X_full))

    results = []
    for name, X_sub in configs:
        runner = AblationRunner(X_sub, y)
        config = AblationConfig(
            name=f"{cancer_type}:{name}",
            encoder_level=0,
            n_estimators=100,
            n_folds=5,
        )
        result = runner.run_config(config)
        results.append(result)

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Per-type ablation study")
    parser.add_argument("--n-types", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 70)
    print("PER-CANCER-TYPE ABLATION: Mechanistic Feature Value")
    print("=" * 70)

    # ── Load MAF ───────────────────────────────────────────────────
    maf_uncompressed = DATA_DIR / "tcga" / "mc3.v0.2.8.PUBLIC.maf"
    maf_compressed = DATA_DIR / "tcga" / "mc3.v0.2.8.PUBLIC.maf.gz"
    maf_path = maf_uncompressed if maf_uncompressed.exists() else maf_compressed
    if not maf_path.exists():
        print("SKIP: MC3 MAF not found")
        return

    # ── Load barcode master ────────────────────────────────────────
    master_path = DATA_DIR / "tcga" / "barcode_master.tsv"
    if not master_path.exists():
        print("SKIP: barcode_master.tsv not found")
        return

    master_df = pd.read_csv(master_path, sep="\t")
    cancer_barcodes = master_df.groupby("cancer_type")["barcode"].apply(set).to_dict()
    type_counts = {ct: len(bcs) for ct, bcs in cancer_barcodes.items()}
    top_types = sorted(type_counts, key=type_counts.get, reverse=True)[:args.n_types]

    # ── Load ground truth ──────────────────────────────────────────
    gt_dir = DATA_DIR / "ground_truth"
    consensus = build_consensus_ground_truth(
        intogen_path=gt_dir / "intogen" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv",
        ncg_path=gt_dir / "NCG_cancerdrivers_annotation.tsv",
        cancermine_path=gt_dir / "cancermine_collated.tsv",
        min_sources=2,
    )
    labels = consensus.as_labels()
    print(f"Ground truth: {len(labels)} consensus drivers")

    # ── Load optional features ─────────────────────────────────────
    depmap_features: dict[str, float] = {}
    depmap_path = DATA_DIR / "depmap" / "CRISPRGeneEffect.csv"
    if depmap_path.exists():
        df_dep = pd.read_csv(depmap_path, index_col=0)
        for col in df_dep.columns:
            gene = col.split(" (")[0] if " (" in col else col
            depmap_features[gene] = float(df_dep[col].mean())
        print(f"DepMap: {len(depmap_features)} genes")
        del df_dep

    gnomad_features: dict[str, float] = {}
    gnomad_path = DATA_DIR / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv"
    if gnomad_path.exists():
        df_gn = pd.read_csv(gnomad_path, sep="\t")
        gene_col = "gene" if "gene" in df_gn.columns else df_gn.columns[0]
        for cand in ("oe_lof_upper", "lof.oe_ci.upper", "oe_lof_upper_bin"):
            if cand in df_gn.columns:
                for _, row in df_gn.iterrows():
                    gene = str(row[gene_col]).strip()
                    val = row[cand]
                    if pd.notna(val):
                        gnomad_features[gene] = float(val)
                break
        print(f"gnomAD: {len(gnomad_features)} genes")
        del df_gn

    # ── Per-type ablation ──────────────────────────────────────────
    all_rows: list[dict] = []

    for cancer_type in top_types:
        barcodes = cancer_barcodes[cancer_type]
        t0 = time.time()

        print(f"\n{'─' * 50}")
        print(f"{cancer_type} ({len(barcodes)} samples)")
        print(f"{'─' * 50}")

        # Parse MAF
        alt_matrix = parse_maf_to_matrix(maf_path, min_samples=2, sample_barcodes=barcodes)
        n_genes = alt_matrix.n_genes
        print(f"  {n_genes} genes x {alt_matrix.n_samples} samples")

        if n_genes < 20:
            print("  SKIP: too few genes")
            continue

        # Build co-alteration graph
        multi_view = build_coalteration_multi_view(
            alt_matrix.data, alt_matrix.genes, alt_matrix.n_samples,
        )

        # Assemble features
        additional = {}
        if depmap_features:
            additional["depmap_essentiality"] = depmap_features
        if gnomad_features:
            additional["gnomad_loeuf"] = gnomad_features

        assembler = FeatureAssembler(
            alteration_matrix=alt_matrix,
            multi_view=multi_view,
            additional_features=additional,
        )
        X, gene_list = assembler.assemble(standardize=True)
        y = np.array([labels.get(g, 0) for g in gene_list])
        omics_idx, topo_idx = assembler.get_bimodal_indices()

        n_pos = int(y.sum())
        print(f"  Features: {assembler.feature_names}")
        print(f"  Positives: {n_pos}/{len(y)}")

        if n_pos < 5:
            print("  SKIP: too few positives for cross-validation")
            del alt_matrix, multi_view, assembler
            gc.collect()
            continue

        # Run ablations
        results = run_ablations_for_type(
            cancer_type, X, y,
            assembler.feature_names, omics_idx, topo_idx,
        )

        # Print & collect
        print(f"\n  {'Config':<18} {'AUPRC':>7} {'AUROC':>7} {'P@10':>6} {'P@50':>6}")
        print(f"  {'-'*46}")

        mut_auprc = None
        for r in results:
            config_name = r.name.split(":", 1)[1]
            if mut_auprc is None:
                mut_auprc = r.auprc
            lift = ((r.auprc / mut_auprc - 1) * 100) if mut_auprc > 0 else 0
            print(f"  {config_name:<18} {r.auprc:>7.4f} {r.auroc:>7.4f} "
                  f"{r.precision_at_10:>6.2f} {r.precision_at_50:>6.2f}  "
                  f"{'(baseline)' if config_name == 'mut_only' else f'+{lift:.0f}%'}")

            all_rows.append({
                "cancer_type": cancer_type,
                "n_genes": n_genes,
                "n_samples": alt_matrix.n_samples,
                "n_positives": n_pos,
                "config": config_name,
                "auprc": f"{r.auprc:.6f}",
                "auroc": f"{r.auroc:.6f}",
                "precision_at_10": f"{r.precision_at_10:.4f}",
                "precision_at_50": f"{r.precision_at_50:.4f}",
                "lift_auprc_pct": f"{lift:.1f}",
            })

        dt = time.time() - t0
        print(f"  Time: {dt:.1f}s")

        del alt_matrix, multi_view, assembler
        gc.collect()

    # ── Save ───────────────────────────────────────────────────────
    df = pd.DataFrame(all_rows)
    df.to_csv(RESULTS_DIR / "per_type_ablation.tsv", sep="\t", index=False)
    print(f"\n{'=' * 70}")
    print(f"Saved: per_type_ablation.tsv ({len(all_rows)} rows)")

    # Summary: types sorted by mutation load (n_genes)
    print(f"\n{'=' * 70}")
    print("MECHANISTIC FEATURE VALUE BY MUTATION LOAD")
    print(f"{'=' * 70}")
    print(f"{'Type':>12} {'Genes':>6} {'mut_only':>10} {'full':>10} {'Lift':>7}")
    print(f"{'-' * 50}")

    for ct in sorted(top_types, key=lambda c: type_counts[c]):
        ct_rows = df[df["cancer_type"] == ct]
        if ct_rows.empty:
            continue
        mut_row = ct_rows[ct_rows["config"] == "mut_only"]
        full_row = ct_rows[ct_rows["config"] == "full"]
        if mut_row.empty or full_row.empty:
            continue
        mut_auprc = float(mut_row.iloc[0]["auprc"])
        full_auprc = float(full_row.iloc[0]["auprc"])
        n_genes = int(mut_row.iloc[0]["n_genes"])
        lift = ((full_auprc / mut_auprc - 1) * 100) if mut_auprc > 0 else 0
        print(f"{ct:>12} {n_genes:>6} {mut_auprc:>10.4f} {full_auprc:>10.4f} {lift:>+6.0f}%")

    print(f"\nTotal time: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
