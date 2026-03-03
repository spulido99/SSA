#!/usr/bin/env python
"""CLI launcher for SageMaker SSA-ME pipeline.

Usage:
    python scripts/cloud/launch_pipeline.py sync          # upload data to S3
    python scripts/cloud/launch_pipeline.py run           # start pipeline
    python scripts/cloud/launch_pipeline.py run --from-step Run_BRCA  # selective execution
    python scripts/cloud/launch_pipeline.py status        # check progress
    python scripts/cloud/launch_pipeline.py results       # download results
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SSA-ME SageMaker Pipeline Launcher")
    sub = parser.add_subparsers(dest="command", required=True)

    # sync
    sync_p = sub.add_parser("sync", help="Upload data to S3")
    sync_p.add_argument("--bucket", required=True, help="S3 bucket name")
    sync_p.add_argument("--prefix", default="ssa-me", help="S3 key prefix")

    # run
    run_p = sub.add_parser("run", help="Start pipeline")
    run_p.add_argument("--bucket", required=True, help="S3 bucket name")
    run_p.add_argument("--prefix", default="ssa-me", help="S3 key prefix")
    run_p.add_argument("--role", required=True, help="SageMaker IAM role ARN")
    run_p.add_argument("--image-uri", required=True, help="Docker image URI")
    run_p.add_argument("--encoder-level", type=int, default=0, choices=[0, 1])
    run_p.add_argument("--n-types", type=int, default=10)
    run_p.add_argument("--n-iter", type=int, default=200)
    run_p.add_argument("--n-ants", type=int, default=20)
    run_p.add_argument("--from-step", default=None, help="Re-run from this step onward")
    run_p.add_argument("--no-spot", action="store_true", help="Disable spot instances")

    # status
    sub.add_parser("status", help="Check pipeline status")

    # results
    res_p = sub.add_parser("results", help="Download results from S3")
    res_p.add_argument("--bucket", required=True, help="S3 bucket name")
    res_p.add_argument("--prefix", default="ssa-me", help="S3 key prefix")
    res_p.add_argument("--output-dir", default="results/cloud", help="Local output dir")

    return parser.parse_args()


def cmd_sync(args: argparse.Namespace) -> None:
    """Upload local data to S3."""
    from scripts.cloud.s3_sync import sync_to_s3

    data_dir = PROJECT_ROOT / "data" / "raw"
    if not data_dir.exists():
        print(f"ERROR: Data directory not found: {data_dir}")
        sys.exit(1)
    print(f"Syncing {data_dir} -> s3://{args.bucket}/{args.prefix}/data/")
    sync_to_s3(data_dir, args.bucket, f"{args.prefix}/data")
    print("Sync complete.")


def cmd_run(args: argparse.Namespace) -> None:
    """Start the SageMaker pipeline."""
    import pandas as pd

    from scripts.cloud.sagemaker_pipeline import build_pipeline_definition, create_and_run_pipeline

    # Load barcode master to determine cancer types and sizes
    master_path = PROJECT_ROOT / "data" / "raw" / "tcga" / "barcode_master.tsv"
    if not master_path.exists():
        print(f"ERROR: Barcode master not found: {master_path}")
        sys.exit(1)

    master_df = pd.read_csv(master_path, sep="\t")
    cancer_barcodes = master_df.groupby("cancer_type")["barcode"].apply(set).to_dict()
    type_counts = {ct: len(bcs) for ct, bcs in cancer_barcodes.items()}
    top_types = sorted(type_counts, key=type_counts.get, reverse=True)[:args.n_types]

    cancer_types = [
        {"name": ct, "n_samples": type_counts[ct], "n_genes": 5000}
        for ct in top_types
    ]

    definition = build_pipeline_definition(
        role=args.role,
        bucket=args.bucket,
        prefix=args.prefix,
        cancer_types=cancer_types,
        encoder_level=args.encoder_level,
        n_iter=args.n_iter,
        n_ants=args.n_ants,
        image_uri=args.image_uri,
        enable_spot=not args.no_spot,
    )

    print(f"Pipeline: {len(definition['steps'])} steps")
    for step in definition["steps"]:
        deps = ", ".join(step["depends_on"]) or "none"
        print(f"  {step['name']}: {step['instance_type']} (depends: {deps})")

    execution_arn = create_and_run_pipeline(
        definition,
        from_step=args.from_step,
    )
    print(f"\nExecution started: {execution_arn}")


def cmd_status(_args: argparse.Namespace) -> None:
    """Check pipeline execution status."""
    try:
        import boto3

        sm = boto3.client("sagemaker")
        response = sm.list_pipeline_executions(PipelineName="ssa-me-pan-cancer")
        executions = response.get("PipelineExecutionSummaries", [])
        if not executions:
            print("No pipeline executions found.")
            return
        for ex in executions[:5]:
            print(f"  {ex['PipelineExecutionArn']}: {ex['PipelineExecutionStatus']}")
    except ImportError:
        print("ERROR: boto3 required. Install with: pip install boto3")
    except Exception as e:
        print(f"ERROR: {e}")


def cmd_results(args: argparse.Namespace) -> None:
    """Download results from S3."""
    try:
        import boto3

        s3 = boto3.client("s3")
        prefix = f"{args.prefix}/final/"
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        response = s3.list_objects_v2(Bucket=args.bucket, Prefix=prefix)
        for obj in response.get("Contents", []):
            key = obj["Key"]
            filename = key.split("/")[-1]
            if filename:
                local_path = output_dir / filename
                print(f"Downloading {key} -> {local_path}")
                s3.download_file(args.bucket, key, str(local_path))
        print(f"Results downloaded to {output_dir}")
    except ImportError:
        print("ERROR: boto3 required. Install with: pip install boto3")
    except Exception as e:
        print(f"ERROR: {e}")


def main() -> None:
    args = parse_args()
    commands = {
        "sync": cmd_sync,
        "run": cmd_run,
        "status": cmd_status,
        "results": cmd_results,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
