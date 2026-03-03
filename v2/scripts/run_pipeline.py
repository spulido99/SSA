#!/usr/bin/env python
"""Run the full SSA-ME v2 pipeline from command line."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np

from ssame.data.alteration_matrix import load_m2
from ssame.network.multi_view_graph import build_multi_view
from ssame.pipeline.pipeline import Pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description="SSA-ME v2 pipeline")
    parser.add_argument(
        "--alteration-matrix", "-m", required=True, help="Path to .m2 alteration file"
    )
    parser.add_argument(
        "--network", "-n", required=True, help="Path to network edge list (TSV)"
    )
    parser.add_argument(
        "--labels", "-l", required=True,
        help="Path to labels file (gene<TAB>1 or gene<TAB>0, one per line)",
    )
    parser.add_argument("--output", "-o", default="results.json", help="Output file")
    parser.add_argument("--n-estimators", type=int, default=100)
    parser.add_argument("--n-ants", type=int, default=50)
    parser.add_argument("--n-iterations", type=int, default=1000)
    parser.add_argument("--subnetwork-size", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-workers", type=int, default=0, help="Thread limit (0=auto)")
    parser.add_argument(
        "--resume", nargs="?", const="", default=None, metavar="STEP",
        help="Resume from checkpoints. Optional step name to fork from "
             "(keeps steps before STEP, re-runs from STEP onward).",
    )
    parser.add_argument(
        "--fresh", action="store_true",
        help="Clear all checkpoints before running.",
    )
    parser.add_argument(
        "--mode", choices=["lite", "standard", "full"], default=None,
        help="Apply parameter preset (overrides --n-ants/--n-iterations/etc.)",
    )
    args = parser.parse_args()

    # Apply preset if given
    from ssame.config import PRESETS
    if args.mode:
        preset = PRESETS[args.mode]
        args.n_ants = preset.get("n_ants", args.n_ants)
        args.n_iterations = preset.get("n_iter", args.n_iterations)
        args.n_estimators = preset.get("n_estimators", args.n_estimators)

    # Load data
    print(f"Loading alteration matrix from {args.alteration_matrix}...")
    alt_matrix = load_m2(Path(args.alteration_matrix))

    print(f"Loading network from {args.network}...")
    graph = nx.read_edgelist(args.network, delimiter="\t")
    multi_view = build_multi_view({"ppi": graph})

    print(f"Loading labels from {args.labels}...")
    labels: dict[str, int] = {}
    with open(args.labels) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                labels[parts[0]] = int(parts[1])

    # Checkpoint setup
    checkpoint_dir = None
    if args.resume is not None:
        checkpoint_dir = Path(args.output).parent / ".checkpoints"
        if args.resume:
            from ssame.utils.checkpoint import CheckpointManager
            ckpt = CheckpointManager(checkpoint_dir)
            ckpt.invalidate_from(args.resume)
            print(f"Forking: invalidated from '{args.resume}' onward")
    if args.fresh and checkpoint_dir and checkpoint_dir.exists():
        import shutil
        shutil.rmtree(checkpoint_dir)
        print("Cleared all checkpoints")

    # Run pipeline
    print("Running pipeline...")
    pipeline = Pipeline(
        alteration_matrix=alt_matrix,
        multi_view=multi_view,
        labels=labels,
        encoder_n_estimators=args.n_estimators,
        search_n_ants=args.n_ants,
        search_n_iterations=args.n_iterations,
        search_subnetwork_size=args.subnetwork_size,
        seed=args.seed,
        checkpoint_dir=checkpoint_dir,
        n_workers=args.n_workers,
    )
    result = pipeline.run()

    # Write output
    output = {
        "ranked_genes": result.ranked_genes[:100],
        "n_modules": len(result.modules),
        "top_modules": [
            {"genes": [g.symbol for g in m.genes], "score": s}
            for m, s in zip(result.modules[:20], result.module_scores[:20], strict=False)
        ],
    }

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Results written to {args.output}")
    print(f"Top 10 ranked genes:")
    for gene, score in result.ranked_genes[:10]:
        print(f"  {gene}: {score:.4f}")


if __name__ == "__main__":
    main()
