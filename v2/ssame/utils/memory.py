"""Memory monitoring utility."""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def log_memory(step_name: str) -> None:
    """Log current memory usage. Warns if < 2GB available."""
    try:
        import psutil

        mem = psutil.virtual_memory()
        logger.info(
            "[MEMORY] %s: %.1fGB used, %.1fGB free",
            step_name,
            mem.used / 1e9,
            mem.available / 1e9,
        )
        if mem.available < 2e9:
            logger.warning("[MEMORY] WARNING: < 2GB available")
    except ImportError:
        pass  # psutil is optional
