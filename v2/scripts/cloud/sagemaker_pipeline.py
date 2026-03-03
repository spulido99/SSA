"""SageMaker Pipeline definition for pan-cancer SSA-ME analysis.

Multi-step pipeline:
  Step 1: ParseMAF — parse MAF once into per-cancer-type alteration matrices
  Step 2a-2j: RunPerCancerType — run pipeline per cancer type (parallel)
  Step 3: Aggregate — Borda aggregation + evaluation

Cancer types run in parallel on separate instances, each right-sized.
"""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)

# Instance sizing thresholds
_GPU_INSTANCE = "ml.g4dn.xlarge"  # 4 vCPU, 16GB, T4 ($0.526/hr)
_CPU_INSTANCES = {
    "large": "ml.m5.large",      # 2 vCPU, 8GB  ($0.115/hr)
    "xlarge": "ml.m5.xlarge",    # 4 vCPU, 16GB ($0.23/hr)
    "2xlarge": "ml.m5.2xlarge",  # 8 vCPU, 32GB ($0.461/hr)
}


def select_instance(n_samples: int, n_genes: int, encoder_level: int) -> str:
    """Select smallest instance that fits the cancer type's data.

    Args:
        n_samples: Number of samples in this cancer type.
        n_genes: Number of genes in the alteration matrix.
        encoder_level: 0 (CPU) or 1 (GPU).

    Returns:
        SageMaker instance type string.
    """
    if encoder_level == 1:
        return _GPU_INSTANCE

    # Estimate memory: gene×sample matrix × 4 bytes × ~3x for intermediates
    mem_estimate_gb = (n_genes * n_samples * 4) / 1e9 * 3
    if mem_estimate_gb > 8:
        return _CPU_INSTANCES["2xlarge"]
    if mem_estimate_gb > 4:
        return _CPU_INSTANCES["xlarge"]
    return _CPU_INSTANCES["large"]


def build_pipeline_definition(
    *,
    role: str,
    bucket: str,
    prefix: str = "ssa-me",
    cancer_types: list[dict[str, Any]],
    encoder_level: int = 0,
    n_iter: int = 200,
    n_ants: int = 20,
    image_uri: str | None = None,
    enable_spot: bool = True,
    cache_days: int = 30,
) -> dict[str, Any]:
    """Build a SageMaker Pipeline definition dict.

    This returns a serializable dict describing the pipeline structure.
    Actual SageMaker SDK objects are created in `create_and_run_pipeline()`.

    Args:
        role: IAM role ARN for SageMaker.
        bucket: S3 bucket for data and outputs.
        prefix: S3 key prefix.
        cancer_types: List of dicts with 'name', 'n_samples', 'n_genes'.
        encoder_level: 0 (XGBoost) or 1 (Graph Transformer).
        n_iter: ACO iterations.
        n_ants: Number of ants.
        image_uri: Docker image URI for processing jobs.
        enable_spot: Use spot instances for per-cancer steps.
        cache_days: Step cache expiry in days.

    Returns:
        Pipeline definition dict with steps, parameters, and instance configs.
    """
    steps = []

    # Step 1: ParseMAF
    steps.append({
        "name": "ParseMAF",
        "instance_type": _CPU_INSTANCES["xlarge"],
        "script": "scripts/cloud/steps/parse_maf.py",
        "inputs": [f"s3://{bucket}/{prefix}/data/tcga/"],
        "outputs": [f"s3://{bucket}/{prefix}/parsed/"],
        "max_runtime_seconds": 3600,
        "depends_on": [],
    })

    # Steps 2a-2j: Per-cancer-type pipeline (all depend on ParseMAF only)
    for ct_info in cancer_types:
        ct_name = ct_info["name"]
        instance = select_instance(
            n_samples=ct_info["n_samples"],
            n_genes=ct_info["n_genes"],
            encoder_level=encoder_level,
        )
        steps.append({
            "name": f"Run_{ct_name}",
            "instance_type": instance,
            "script": "scripts/cloud/steps/run_cancer_type.py",
            "inputs": [f"s3://{bucket}/{prefix}/parsed/{ct_name}/"],
            "outputs": [f"s3://{bucket}/{prefix}/results/{ct_name}/"],
            "max_runtime_seconds": 14400,
            "depends_on": ["ParseMAF"],
            "use_spot": enable_spot,
            "env": {
                "CANCER_TYPE": ct_name,
                "ENCODER_LEVEL": str(encoder_level),
                "N_ITER": str(n_iter),
                "N_ANTS": str(n_ants),
            },
        })

    # Step 3: Aggregate
    agg_depends = [f"Run_{ct['name']}" for ct in cancer_types]
    steps.append({
        "name": "Aggregate",
        "instance_type": _CPU_INSTANCES["large"],
        "script": "scripts/cloud/steps/aggregate.py",
        "inputs": [f"s3://{bucket}/{prefix}/results/"],
        "outputs": [f"s3://{bucket}/{prefix}/final/"],
        "max_runtime_seconds": 3600,
        "depends_on": agg_depends,
    })

    return {
        "pipeline_name": "ssa-me-pan-cancer",
        "role": role,
        "bucket": bucket,
        "prefix": prefix,
        "image_uri": image_uri,
        "cache_days": cache_days,
        "parameters": {
            "EncoderLevel": encoder_level,
            "NumCancerTypes": len(cancer_types),
            "NumACOIterations": n_iter,
            "NumAnts": n_ants,
        },
        "steps": steps,
    }


def create_and_run_pipeline(
    definition: dict[str, Any],
    *,
    from_step: str | None = None,
) -> str:
    """Create and execute a SageMaker Pipeline from definition.

    Args:
        definition: Pipeline definition from build_pipeline_definition().
        from_step: If set, selectively execute from this step onward
                   (reuses cached outputs of earlier steps).

    Returns:
        Execution ARN string.
    """
    try:
        import sagemaker
        from sagemaker.processing import ProcessingInput, ProcessingOutput, ScriptProcessor
        from sagemaker.workflow.pipeline import Pipeline
        from sagemaker.workflow.steps import CacheConfig, ProcessingStep
    except ImportError as e:
        msg = "sagemaker SDK required. Install with: pip install sagemaker"
        raise ImportError(msg) from e

    role = definition["role"]
    image_uri = definition["image_uri"]
    cache_config = CacheConfig(
        enable_caching=True,
        expire_after=f"P{definition['cache_days']}D",
    )

    sm_steps = []
    step_map: dict[str, ProcessingStep] = {}

    for step_def in definition["steps"]:
        processor = ScriptProcessor(
            role=role,
            image_uri=image_uri,
            instance_type=step_def["instance_type"],
            instance_count=1,
            command=["python3"],
            max_runtime_in_seconds=step_def["max_runtime_seconds"],
            env=step_def.get("env", {}),
        )

        inputs = [
            ProcessingInput(source=src, destination="/opt/ml/processing/input/")
            for src in step_def["inputs"]
        ]
        outputs = [
            ProcessingOutput(source="/opt/ml/processing/output/", destination=dst)
            for dst in step_def["outputs"]
        ]

        sm_step = ProcessingStep(
            name=step_def["name"],
            processor=processor,
            code=step_def["script"],
            inputs=inputs,
            outputs=outputs,
            cache_config=cache_config,
        )

        # Wire dependencies
        for dep_name in step_def["depends_on"]:
            if dep_name in step_map:
                sm_step.add_depends_on([step_map[dep_name]])

        step_map[step_def["name"]] = sm_step
        sm_steps.append(sm_step)

    pipeline = Pipeline(
        name=definition["pipeline_name"],
        steps=sm_steps,
        sagemaker_session=sagemaker.Session(),
    )

    pipeline.upsert(role_arn=role)
    logger.info("Pipeline created/updated: %s", definition["pipeline_name"])

    if from_step:
        # Selective execution: reuse earlier step outputs
        from sagemaker.workflow.selective_execution_config import SelectiveExecutionConfig

        selected = [s["name"] for s in definition["steps"]
                    if s["name"] == from_step
                    or any(from_step in dep for dep in [s.get("depends_on", [])])]
        # Include the from_step and all downstream
        from_idx = next(
            i for i, s in enumerate(definition["steps"]) if s["name"] == from_step
        )
        selected = [s["name"] for s in definition["steps"][from_idx:]]

        execution = pipeline.start(
            selective_execution_config=SelectiveExecutionConfig(
                selected_steps=selected,
                reference_latest_execution=True,
            )
        )
    else:
        execution = pipeline.start()

    logger.info("Pipeline execution started: %s", execution.arn)
    return execution.arn
