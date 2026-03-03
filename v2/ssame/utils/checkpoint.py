"""Checkpoint manager for pipeline resume support."""

from __future__ import annotations

import logging
import shutil
import tempfile
from pathlib import Path
from typing import Any

import joblib

logger = logging.getLogger(__name__)

# Ordered pipeline steps for invalidate_from support.
STEP_ORDER = ("features", "scores", "modules", "result")


class CheckpointManager:
    """Save and load intermediate pipeline results for resume support.

    Uses joblib for efficient numpy/scipy serialization.
    Writes atomically (temp file + rename) to prevent corruption.
    """

    def __init__(self, checkpoint_dir: Path, enabled: bool = True) -> None:
        self._dir = checkpoint_dir
        self._enabled = enabled
        if self._enabled:
            self._dir.mkdir(parents=True, exist_ok=True)

    def save(self, key: str, data: Any) -> None:
        """Save data to checkpoint. Atomic write via temp+rename."""
        if not self._enabled:
            return
        self._dir.mkdir(parents=True, exist_ok=True)
        target = self._dir / f"{key}.pkl"
        fd, tmp_path = tempfile.mkstemp(suffix=".tmp", dir=str(self._dir))
        try:
            import os

            os.close(fd)
            joblib.dump(data, tmp_path)
            Path(tmp_path).replace(target)
        except BaseException:
            Path(tmp_path).unlink(missing_ok=True)
            raise
        size_mb = target.stat().st_size / 1e6
        logger.info("[CHECKPOINT] Saved %s (%.1fMB) -> %s", key, size_mb, target)

    def load(self, key: str) -> Any | None:
        """Load checkpoint data. Returns None if missing or disabled."""
        if not self._enabled:
            return None
        path = self._dir / f"{key}.pkl"
        if not path.exists():
            return None
        data = joblib.load(path)
        logger.info("[CHECKPOINT] Loaded %s from %s", key, path)
        return data

    def exists(self, key: str) -> bool:
        """Check if a checkpoint exists."""
        if not self._enabled:
            return False
        return (self._dir / f"{key}.pkl").exists()

    def clear(self) -> None:
        """Delete all checkpoints."""
        if self._dir.exists():
            shutil.rmtree(self._dir)
            logger.info("[CHECKPOINT] Cleared all checkpoints in %s", self._dir)

    def invalidate_from(self, step: str) -> None:
        """Delete the given step and all downstream checkpoints.

        Step order: features -> scores -> modules -> result.
        Unknown keys are deleted individually.
        """
        if step in STEP_ORDER:
            idx = STEP_ORDER.index(step)
            for s in STEP_ORDER[idx:]:
                path = self._dir / f"{s}.pkl"
                if path.exists():
                    path.unlink()
                    logger.info("[CHECKPOINT] Invalidated %s", s)
        else:
            # For non-standard keys (e.g., cancer_BRCA), delete just that one
            path = self._dir / f"{step}.pkl"
            if path.exists():
                path.unlink()
                logger.info("[CHECKPOINT] Invalidated %s", step)
