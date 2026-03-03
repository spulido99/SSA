"""Tests for CheckpointManager."""

import numpy as np
import pytest

from ssame.utils.checkpoint import CheckpointManager


@pytest.fixture
def ckpt_dir(tmp_path):
    """Return a temporary checkpoint directory."""
    return tmp_path / ".checkpoints"


@pytest.fixture
def manager(ckpt_dir):
    """Return an enabled CheckpointManager."""
    return CheckpointManager(checkpoint_dir=ckpt_dir)


class TestSaveAndLoad:
    def test_save_and_load_numpy_array(self, manager):
        arr = np.array([1.0, 2.0, 3.0])
        manager.save("features", arr)
        loaded = manager.load("features")
        np.testing.assert_array_equal(loaded, arr)

    def test_save_and_load_dict(self, manager):
        data = {"genes": ["TP53", "KRAS"], "scores": [0.9, 0.8]}
        manager.save("scores", data)
        loaded = manager.load("scores")
        assert loaded == data

    def test_save_and_load_tuple(self, manager):
        arr = np.ones((10, 5))
        genes = ["A", "B", "C"]
        manager.save("features", (arr, genes))
        loaded_arr, loaded_genes = manager.load("features")
        np.testing.assert_array_equal(loaded_arr, arr)
        assert loaded_genes == genes


class TestLoadMissing:
    def test_load_missing_returns_none(self, manager):
        assert manager.load("nonexistent") is None


class TestExists:
    def test_exists_false_before_save(self, manager):
        assert not manager.exists("features")

    def test_exists_true_after_save(self, manager):
        manager.save("features", [1, 2, 3])
        assert manager.exists("features")


class TestClear:
    def test_clear_removes_all(self, manager):
        manager.save("features", [1])
        manager.save("scores", [2])
        manager.clear()
        assert not manager.exists("features")
        assert not manager.exists("scores")


class TestDisabled:
    def test_disabled_manager_skips_save(self, ckpt_dir):
        mgr = CheckpointManager(checkpoint_dir=ckpt_dir, enabled=False)
        mgr.save("features", [1, 2, 3])
        assert not mgr.exists("features")

    def test_disabled_manager_load_returns_none(self, ckpt_dir):
        mgr = CheckpointManager(checkpoint_dir=ckpt_dir, enabled=False)
        assert mgr.load("features") is None


class TestAtomicWrite:
    def test_atomic_write_no_partial_file_on_disk(self, manager, ckpt_dir):
        """Save should not leave partial files — uses temp+rename."""
        manager.save("data", np.zeros(1000))
        # After save completes, only the final .pkl should exist
        pkl_files = list(ckpt_dir.glob("*.pkl"))
        assert len(pkl_files) == 1
        assert pkl_files[0].name == "data.pkl"
        # No temp files left behind
        tmp_files = list(ckpt_dir.glob("*.tmp"))
        assert len(tmp_files) == 0


class TestInvalidateFrom:
    def test_invalidate_from_deletes_step_and_downstream(self, manager):
        manager.save("features", [1])
        manager.save("scores", [2])
        manager.save("modules", [3])
        manager.save("result", [4])

        manager.invalidate_from("scores")

        assert manager.exists("features")
        assert not manager.exists("scores")
        assert not manager.exists("modules")
        assert not manager.exists("result")

    def test_invalidate_from_first_step_clears_all(self, manager):
        manager.save("features", [1])
        manager.save("scores", [2])

        manager.invalidate_from("features")

        assert not manager.exists("features")
        assert not manager.exists("scores")

    def test_invalidate_from_last_step_only_deletes_it(self, manager):
        manager.save("features", [1])
        manager.save("scores", [2])
        manager.save("modules", [3])
        manager.save("result", [4])

        manager.invalidate_from("result")

        assert manager.exists("features")
        assert manager.exists("scores")
        assert manager.exists("modules")
        assert not manager.exists("result")


class TestLogging:
    def test_save_logs_message(self, manager, caplog):
        import logging

        with caplog.at_level(logging.INFO):
            manager.save("features", np.zeros(100))

        assert any("[CHECKPOINT]" in r.message for r in caplog.records)
        assert any("features" in r.message for r in caplog.records)
