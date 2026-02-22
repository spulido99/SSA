"""Save pipeline and ablation results to CSV/TSV for paper figures.

Runs after smoke_test.py and run_ablation.py to persist structured results.

Usage:
    python scripts/save_results.py
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

from ssame.data.alteration_matrix import parse_maf_to_matrix
from ssame.data.feature_assembly import FeatureAssembler
from ssame.data.open_ground_truth import build_consensus_ground_truth
from ssame.evaluation.ablation import AblationConfig, AblationRunner
from ssame.network.multi_view_graph import build_multi_view
from ssame.network.pioneer_oncoppi import load_pioneer
from ssame.pipeline.pipeline import Pipeline

RESULTS_DIR = PROJECT_ROOT / "results"


def save_ground_truth_stats(consensus, out_dir: Path) -> None:
    """Save ground truth source overlap statistics."""
    rows = []
    for gene in sorted(consensus.source_count.keys()):
        rows.append({
            "gene": gene,
            "n_sources": consensus.source_count[gene],
            "consensus_role": consensus.consensus_role.get(gene, "unknown"),
            "is_consensus": gene in consensus.positives,
            "is_high_confidence": gene in consensus.high_confidence,
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "ground_truth_genes.tsv", sep="\t", index=False)

    # Summary
    summary = {
        "total_genes_any_source": len(consensus.source_count),
        "consensus_drivers": len(consensus.positives),
        "high_confidence_drivers": len(consensus.high_confidence),
    }
    for source_name, result in consensus.per_source.items():
        summary[f"{source_name}_count"] = len(result.positives)

    with open(out_dir / "ground_truth_summary.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        for k, v in summary.items():
            writer.writerow([k, v])


def save_pipeline_results(result, labels, out_dir: Path) -> None:
    """Save ranked genes and module info."""
    # Ranked genes
    rows = []
    for rank, (gene, score) in enumerate(result.ranked_genes, 1):
        rows.append({
            "rank": rank,
            "gene": gene,
            "score": f"{score:.6f}",
            "is_known_driver": gene in labels,
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "ranked_genes.tsv", sep="\t", index=False)

    # Module summary
    rows = []
    for i, (module, score) in enumerate(zip(result.modules, result.module_scores)):
        gene_symbols = sorted(g.symbol for g in module.genes)
        n_drivers = sum(1 for g in gene_symbols if g in labels)
        rows.append({
            "module_id": i,
            "score": f"{score:.6f}",
            "size": len(gene_symbols),
            "n_known_drivers": n_drivers,
            "genes": ",".join(gene_symbols),
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "modules.tsv", sep="\t", index=False)


def save_ablation_results(results, out_dir: Path) -> None:
    """Save ablation table."""
    rows = []
    for r in results:
        rows.append({
            "config": r.name,
            "auprc": f"{r.auprc:.4f}",
            "auroc": f"{r.auroc:.4f}",
            "precision_at_10": f"{r.precision_at_10:.4f}",
            "precision_at_50": f"{r.precision_at_50:.4f}",
            "recall_at_50": f"{r.recall_at_50:.4f}",
            "per_fold_auprc": ",".join(f"{x:.4f}" for x in r.per_fold_auprc),
            "per_fold_auroc": ",".join(f"{x:.4f}" for x in r.per_fold_auroc),
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "ablation_results.tsv", sep="\t", index=False)


def main() -> None:
    data_dir = PROJECT_ROOT / "data" / "raw"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 60)
    print("SAVING RESULTS TO", RESULTS_DIR)
    print("=" * 60)

    # ── Load data ─────────────────────────────────────────────────
    maf_path = data_dir / "tcga" / "mc3.v0.2.8.PUBLIC.maf.gz"
    if not maf_path.exists():
        print(f"SKIP: MC3 MAF not found at {maf_path}")
        return

    print("\nLoading mutation data (200 samples)...")
    first_chunk = pd.read_csv(
        maf_path, sep="\t", comment="#",
        usecols=["Tumor_Sample_Barcode"], nrows=400_000,
    )
    subset = set(first_chunk["Tumor_Sample_Barcode"].unique()[:200])
    alt_matrix = parse_maf_to_matrix(maf_path, min_samples=2, sample_barcodes=subset)
    print(f"  {alt_matrix.data.shape[0]} genes x {alt_matrix.data.shape[1]} samples")

    print("Loading PIONEER...")
    pioneer_path = data_dir / "networks" / "pioneer" / "pioneer_very_high_human.txt"
    pioneer_graph = load_pioneer(pioneer_path)
    multi_view = build_multi_view({"pioneer": pioneer_graph})

    print("Building ground truth...")
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

    # ── Save ground truth stats ───────────────────────────────────
    print("\nSaving ground truth stats...")
    save_ground_truth_stats(consensus, RESULTS_DIR)

    # ── Run pipeline ──────────────────────────────────────────────
    print("Running pipeline (L0, 200 samples)...")
    pipeline = Pipeline(
        alteration_matrix=alt_matrix,
        multi_view=multi_view,
        labels=labels,
        encoder_level=0,
        encoder_n_estimators=100,
        search_n_ants=10,
        search_n_iterations=50,
        search_subnetwork_size=3,
        seed=42,
    )
    result = pipeline.run()
    save_pipeline_results(result, labels, RESULTS_DIR)
    print(f"  Saved: ranked_genes.tsv ({len(result.ranked_genes)} genes)")
    print(f"  Saved: modules.tsv ({len(result.modules)} modules)")

    # ── Run ablation ──────────────────────────────────────────────
    print("Running ablation study...")

    # Load additional features
    additional = {}
    depmap_path = data_dir / "depmap" / "CRISPRGeneEffect.csv"
    if depmap_path.exists():
        df_dep = pd.read_csv(depmap_path, index_col=0)
        depmap_feat = {}
        for col in df_dep.columns:
            gene = col.split(" (")[0] if " (" in col else col
            depmap_feat[gene] = float(df_dep[col].mean())
        additional["depmap_essentiality"] = depmap_feat

    gnomad_path = data_dir / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv"
    if gnomad_path.exists():
        df_gn = pd.read_csv(gnomad_path, sep="\t")
        gnomad_feat = {}
        for _, row in df_gn.iterrows():
            gene = str(row["gene"]).strip()
            val = row.get("lof.oe_ci.upper")
            if pd.notna(val):
                gnomad_feat[gene] = float(val)
        additional["gnomad_loeuf"] = gnomad_feat

    assembler = FeatureAssembler(
        alteration_matrix=alt_matrix,
        multi_view=multi_view,
        additional_features=additional,
    )
    X, gene_list = assembler.assemble(standardize=True)
    y = np.array([labels.get(g, 0) for g in gene_list])
    feat_names = assembler.feature_names
    _, topo_idx = assembler.get_bimodal_indices()

    mut_cols = [feat_names.index("mutation_frequency"), feat_names.index("log_mutation_count")]

    from ssame.evaluation.ablation import AblationResult
    ablation_results: list[AblationResult] = []

    configs = [
        ("mut_freq_only", X[:, mut_cols]),
        ("+ topology", X[:, mut_cols + topo_idx]),
    ]
    if "depmap_essentiality" in feat_names:
        dc = feat_names.index("depmap_essentiality")
        configs.append(("+ depmap", X[:, mut_cols + [dc]]))
        configs.append(("+ depmap + topo", X[:, mut_cols + [dc] + topo_idx]))
    if "gnomad_loeuf" in feat_names:
        gc = feat_names.index("gnomad_loeuf")
        configs.append(("+ gnomad", X[:, mut_cols + [gc]]))
        configs.append(("+ gnomad + topo", X[:, mut_cols + [gc] + topo_idx]))
    configs.append(("full_model", X))

    for name, X_sub in configs:
        runner = AblationRunner(X_sub, y)
        cfg = AblationConfig(name=name, encoder_level=0, n_estimators=100, n_folds=5)
        ablation_results.append(runner.run_config(cfg))
        print(f"  {name}: AUPRC={ablation_results[-1].auprc:.4f}")

    for n_est in [50, 200]:
        runner = AblationRunner(X, y)
        cfg = AblationConfig(name=f"full_n_est={n_est}", encoder_level=0, n_estimators=n_est, n_folds=5)
        ablation_results.append(runner.run_config(cfg))
        print(f"  full_n_est={n_est}: AUPRC={ablation_results[-1].auprc:.4f}")

    save_ablation_results(ablation_results, RESULTS_DIR)
    print(f"  Saved: ablation_results.tsv ({len(ablation_results)} configs)")

    t_total = time.time() - t_start
    print(f"\nAll results saved to {RESULTS_DIR} in {t_total:.1f}s")


if __name__ == "__main__":
    main()
