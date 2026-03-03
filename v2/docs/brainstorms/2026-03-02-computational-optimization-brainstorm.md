---
date: 2026-03-02
topic: computational-optimization
---

# Computational Optimization: Checkpoints + Performance + AWS

## What We're Building

A 4-phase incremental improvement to make the SSA-ME v2 pipeline runnable on
a laptop with limited RAM (~16GB) and an RTX 3060 (6GB VRAM), with optional
cloud offload for production runs.

**Primary constraint: MEMORY**, not speed. We prefer slower but memory-safe
approaches over faster ones that risk OOM.

## Phase 1: Checkpoint/Resume System

### Design

Add `CheckpointManager` class that saves/loads intermediate results to disk.

**Two levels:**

1. **Pan-cancer level** (`PanCancerRunner`):
   - After each cancer type completes, save its `PipelineResult` to
     `results/.checkpoints/pan_cancer/{cancer_type}.pkl`
   - On resume, skip already-completed types
   - Critical: `gc.collect()` after saving each type to free memory

2. **Pipeline level** (`Pipeline`):
   - Save state after each major step:
     - Step 1: features (X, gene_list) → `features.pkl`
     - Step 3: encoder scores → `scores.pkl`
     - Step 5: ACO modules → `modules.pkl`
   - On resume, load from last completed checkpoint

**Interface:**

```python
class CheckpointManager:
    def __init__(self, checkpoint_dir: Path, enabled: bool = True):
        ...

    def save(self, key: str, data: Any) -> None:
        """Save data to checkpoint file."""

    def load(self, key: str) -> Any | None:
        """Load checkpoint, returns None if not found."""

    def exists(self, key: str) -> bool:
        """Check if checkpoint exists."""

    def clear(self) -> None:
        """Remove all checkpoints (--fresh flag)."""
```

**CLI flags:** `--resume` (default off), `--fresh` (clear checkpoints)

### Memory benefit

By checkpointing after each cancer type and calling `gc.collect()`, the
pan-cancer runner keeps only 1 cancer type in memory at a time instead of
accumulating all 10. This alone could save ~500MB-1GB.

## Phase 2: Memory-First Local Optimizations

### 2a. Streaming MAF with explicit memory release

**Keep current pattern**: read MAF once per cancer type (10 reads total).
This is intentional — it keeps peak memory at ~400MB (one cancer type's
data) instead of ~4GB (all cancer types cached).

**Improvement**: After building the AlterationMatrix for each type, ensure
PyArrow/pandas intermediate objects are explicitly deleted:

```python
matrix = parse_maf_to_matrix(maf_path, sample_barcodes=barcodes)
# PyArrow table is already GC'd inside the function
gc.collect()  # force collection of any remaining fragments
```

### 2b. Sparse co-alteration matrix

`PanCancerRunner._build_coalteration_graph()` computes `data @ data.T`
which creates a dense N×N float32 matrix. For 5000 genes = 100MB.

**Replace with scipy.sparse:**
- Store alteration matrix as `csr_matrix` (int8)
- Compute co-alteration as sparse `(data @ data.T)` → stays sparse
- Build adjacency as sparse bool → pass to MultiViewGraph

**Memory saving**: ~100MB → ~5MB per cancer type (co-alteration graphs
are typically >95% sparse).

### 2c. Chunked BiologicalPriorLoss

`BiologicalPriorLoss.__init__()` computes `alt @ alt.T` eagerly (N×N).
And `.forward()` computes `scores.unsqueeze(0) * scores.unsqueeze(1)` (N×N).

**For N=5000**: each matrix = 100MB. Two = 200MB during forward pass.
**For N=15000** (full TCGA): each = 900MB. Unworkable on laptop.

**Solution**: Chunk the loss computation:
- Don't store full `alt @ alt.T`; compute blocks on-the-fly
- In `.forward()`, loop over gene chunks of size 1000
- Trades ~5x slower L1 training for ~10x less memory

### 2d. Fix O(n²) gene index lookup

`pipeline.py:161` uses `self.alteration_matrix.genes.index(gene)` which
is O(n) per call, O(n²) total. Replace with a pre-built dict.

**Memory impact**: negligible (small dict), but prevents O(n²) CPU burn.

### 2e. NumPy-vectorized pheromone operations

Convert `_pheromone` from `dict[str, float]` to a NumPy array indexed by
node position. Vectorize `_evaporate()` and `_update_pheromone()`.

**Memory impact**: neutral (array vs dict uses similar memory).
**Speed impact**: ~3-5x faster ACO.

### 2f. GPU for L1 encoder (memory-aware)

Add device selection with VRAM guard:
- Auto-detect CUDA availability
- If VRAM < required (estimate from tensor sizes), stay on CPU
- Use `torch.cuda.empty_cache()` aggressively between steps
- Only move encoder + current batch to GPU, not all data

For typical cancer type (3000 genes, hidden=64):
- Model: ~200KB
- Features: ~1MB
- BiologicalPriors: needs chunking (see 2c)
- **Total VRAM: ~50MB** — well within 6GB

### 2g. Early stopping for encoder

Stop L1 training if loss doesn't improve for 10 epochs. Saves 30-50% of
the default 100 epochs.

### 2h. Replace iterrows() with vectorized pandas

In `run_pan_cancer.py:84` and gnomAD loaders. Use `groupby()` or
`to_dict()` instead of `iterrows()`.

### 2i. Memory presets

Add `--mode` flag:
- `lite`: n_iter=50, n_ants=10, n_estimators=50 (fast dev/debug)
- `standard`: n_iter=200, n_ants=20, n_estimators=100 (laptop production)
- `full`: n_iter=1000, n_ants=50, n_estimators=200 (cloud/server)

## Phase 3: Serial/Low-Memory Mode

### 3a. Explicit garbage collection between pipeline steps

After each major step in `Pipeline.run()`, delete intermediate objects
and call `gc.collect()`. If GPU active, call `torch.cuda.empty_cache()`.

### 3b. Thread/process limiting

Set `torch.set_num_threads(2)` and limit XGBoost `nthread=2` to avoid
CPU thrashing on laptop. Configurable via `--n-workers` flag.

### 3c. Memory monitoring

Add optional memory logging (via `psutil`) that warns when available
memory drops below 2GB:

```
[MEMORY] Step 3/6: 4.2GB used, 3.1GB available
[MEMORY] WARNING: Available memory below 2GB threshold
```

## Phase 4: AWS SageMaker (Optional)

### Architecture

```
Local laptop                          AWS
┌────────────────┐                   ┌──────────────────┐
│ Prepare data   │ ─── S3 upload ──► │ SageMaker        │
│ (MAF, labels,  │                   │ Processing Job   │
│  networks)     │                   │                  │
│                │ ◄── S3 download ─ │ Run pipeline     │
│ Analyze results│                   │ Save checkpoints │
└────────────────┘                   └──────────────────┘
```

### Implementation

1. **`scripts/sagemaker_sync.py`**: smart sync local data → S3
   - Computes MD5/ETag hash of local files before uploading
   - Compares with S3 object ETags — skips unchanged files
   - The 3GB MAF only uploads once, persists in S3 indefinitely
   - Small config/label files re-upload only if modified
   - Uses `boto3` multipart upload for large files (>100MB)
   - Prints: `[S3] mc3.maf: already synced (ETag match), skipping`
2. **`scripts/sagemaker_run.py`**: launches Processing Job
3. **`scripts/sagemaker_results.py`**: downloads results from S3

### Instance recommendations (us-east-1)

| Use case | Instance | vCPU | RAM | GPU | $/hr |
|----------|----------|------|-----|-----|------|
| Pan-cancer L0 | ml.m5.xlarge | 4 | 16GB | - | $0.23 |
| Pan-cancer L1 | ml.g4dn.xlarge | 4 | 16GB | T4 16GB | $0.526 |
| Full pipeline | ml.m5.2xlarge | 8 | 32GB | - | $0.461 |

### Cost estimates

| Workload | Instance | Est. time | Est. cost |
|----------|----------|-----------|-----------|
| Pan-cancer 10 types L0 | ml.m5.xlarge | 2-3h | $0.50-0.70 |
| Pan-cancer 10 types L1 | ml.g4dn.xlarge | 1-2h | $0.53-1.05 |
| Full ablation study | ml.m5.xlarge | 1h | $0.23 |
| **Total full run** | | | **$1.30-2.00** |

### S3 data management

- **Bucket**: `s3://ssa-me-v2-data/` (private, versioned)
- **Layout**:
  ```
  s3://ssa-me-v2-data/
    input/
      tcga/mc3.v0.2.8.PUBLIC.maf          # 3GB, uploaded once
      tcga/barcode_master.tsv             # ~1MB
      ground_truth/...                    # ~10MB total
      networks/...                        # ~50MB total
    output/
      checkpoints/                        # per-job checkpoints
      results/                            # final results
  ```
- **Upload policy**: hash-based dedup. File only re-uploads if local
  MD5 != S3 ETag. The 3GB MAF stays in S3 permanently after first upload.
- **Lifecycle rule**: auto-delete output/checkpoints/ after 30 days
  to control storage costs.

### Cost controls

- Use Spot Instances (60-80% discount, with checkpoint/resume to handle interruptions)
- Set `MaxRuntimeInSeconds` to prevent runaway costs
- S3 storage: ~$0.07/month for 3GB data tier
- No persistent endpoints (Processing Jobs only)
- Smart sync prevents re-uploading unchanged data

## Why This Approach

Memory-first design because:
1. Laptop RAM is the binding constraint (not CPU, not GPU, not disk)
2. Checkpoints solve both crash-recovery AND memory (process one type at a time)
3. GPU optimization is secondary — L0 (XGBoost) doesn't use GPU at all
4. AWS is insurance, not primary strategy — optimize locally first

## Key Decisions

- **Keep 10 MAF reads**: intentional memory trade. 10 × 30s I/O is better
  than 1 × 4GB in memory
- **Sparse matrices**: biggest single memory win for co-alteration
- **Chunk BiologicalPriorLoss**: only matters for L1 encoder, can defer
- **Checkpoints before optimization**: safety net enables bolder changes

## Implementation Order

1. CheckpointManager + gc.collect() in pan-cancer loop
2. Sparse co-alteration + gene index dict fix
3. Memory presets (--mode lite/standard/full)
4. NumPy-vectorized ACO pheromone
5. GPU support for L1 (with VRAM guard)
6. iterrows() fixes
7. Memory monitoring
8. AWS SageMaker scripts

## Open Questions

- Exact RAM on laptop? (assumed 16GB)
- Should checkpoints use pickle or joblib? (joblib better for numpy arrays)
- Want AWS Spot Instances? (cheaper but can be interrupted — checkpoints
  make this safe)

## Next Steps

→ Plan implementation with TDD approach per CLAUDE.md requirements
