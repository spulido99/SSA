"""Tests for S3 smart sync (hash-based upload dedup)."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# Add project root so `scripts.cloud.s3_sync` is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from scripts.cloud.s3_sync import compute_md5, sync_to_s3


class TestComputeMd5:
    def test_deterministic(self, tmp_path: Path) -> None:
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        assert compute_md5(f) == compute_md5(f)

    def test_different_content_different_hash(self, tmp_path: Path) -> None:
        f1 = tmp_path / "a.txt"
        f1.write_text("hello")
        f2 = tmp_path / "b.txt"
        f2.write_text("world")
        assert compute_md5(f1) != compute_md5(f2)


class TestS3Sync:
    def test_skips_unchanged_files(self, tmp_path: Path) -> None:
        """Files with matching MD5/ETag should be skipped."""
        (tmp_path / "data.tsv").write_text("hello")
        local_md5 = compute_md5(tmp_path / "data.tsv")

        mock_s3 = MagicMock()
        mock_s3.head_object.return_value = {"ETag": f'"{local_md5}"'}

        sync_to_s3(tmp_path, "my-bucket", "prefix", s3_client=mock_s3)

        mock_s3.upload_file.assert_not_called()

    def test_uploads_new_files(self, tmp_path: Path) -> None:
        """Files not on S3 should be uploaded."""
        (tmp_path / "data.tsv").write_text("hello")

        mock_s3 = MagicMock()
        mock_s3.head_object.side_effect = Exception("not found")

        sync_to_s3(tmp_path, "my-bucket", "prefix", s3_client=mock_s3)

        mock_s3.upload_file.assert_called_once()
        call_args = mock_s3.upload_file.call_args
        assert call_args[0][1] == "my-bucket"
        assert "prefix/data.tsv" in call_args[0][2]

    def test_uploads_changed_files(self, tmp_path: Path) -> None:
        """Files with different MD5 should be re-uploaded."""
        (tmp_path / "data.tsv").write_text("new content")

        mock_s3 = MagicMock()
        mock_s3.head_object.return_value = {"ETag": '"old_hash_mismatch"'}

        sync_to_s3(tmp_path, "my-bucket", "prefix", s3_client=mock_s3)

        mock_s3.upload_file.assert_called_once()

    def test_skips_directories(self, tmp_path: Path) -> None:
        """Directories should be skipped, only files uploaded."""
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        (subdir / "file.txt").write_text("data")

        mock_s3 = MagicMock()
        mock_s3.head_object.side_effect = Exception("not found")

        sync_to_s3(tmp_path, "my-bucket", "prefix", s3_client=mock_s3)

        # Only the file should be uploaded, not the directory
        assert mock_s3.upload_file.call_count == 1
