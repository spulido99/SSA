"""S3 smart sync: hash-based upload deduplication."""

from __future__ import annotations

import hashlib
from pathlib import Path


def compute_md5(path: Path) -> str:
    """Compute MD5 hash of a local file (matches S3 ETag for single-part uploads)."""
    h = hashlib.md5()  # noqa: S324
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def sync_to_s3(
    local_dir: Path,
    bucket: str,
    prefix: str,
    *,
    s3_client: object | None = None,
) -> None:
    """Upload only changed files by comparing MD5 with S3 ETag.

    Args:
        local_dir: Local directory to sync.
        bucket: S3 bucket name.
        prefix: S3 key prefix.
        s3_client: Optional pre-configured boto3 S3 client (for testing).
    """
    if s3_client is None:
        import boto3

        s3_client = boto3.client("s3")
    s3 = s3_client

    for local_file in sorted(local_dir.rglob("*")):
        if local_file.is_dir():
            continue
        s3_key = f"{prefix}/{local_file.relative_to(local_dir)}"
        local_md5 = compute_md5(local_file)
        try:
            response = s3.head_object(Bucket=bucket, Key=s3_key)
            remote_etag = response["ETag"].strip('"')
            if local_md5 == remote_etag:
                print(f"[S3] {s3_key}: unchanged, skipping")
                continue
        except Exception:  # noqa: BLE001
            pass  # file doesn't exist on S3 yet
        size_mb = local_file.stat().st_size / 1e6
        print(f"[S3] Uploading {local_file.name} ({size_mb:.1f}MB)")
        s3.upload_file(str(local_file), bucket, s3_key)
