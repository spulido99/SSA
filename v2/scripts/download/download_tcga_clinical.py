#!/usr/bin/env python
"""Download TCGA barcode-to-cancer-type mapping via GDC API.

Queries the GDC API to map TCGA sample barcodes (e.g. TCGA-A1-0001-01A-...)
to their project/cancer type (e.g. TCGA-BRCA, TCGA-LUAD).

Usage:
    python scripts/download/download_tcga_clinical.py
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "tcga"
OUTPUT_FILE = OUTPUT_DIR / "tcga_barcode_cancer_map.tsv"

GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"
BATCH_SIZE = 5000
MAX_CASES = 25000


def fetch_gdc_cases(offset: int = 0, size: int = BATCH_SIZE) -> dict:
    """Fetch a batch of TCGA cases from GDC API."""
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "project.program.name", "value": "TCGA"}},
        ],
    }
    params = {
        "filters": json.dumps(filters, separators=(",", ":")),
        "fields": "submitter_id,project.project_id",
        "size": str(size),
        "from": str(offset),
        "format": "JSON",
    }

    query_string = urllib.parse.urlencode(params)
    url = f"{GDC_CASES_URL}?{query_string}"

    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.loads(response.read().decode())


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Delete empty file from previous failed run
    if OUTPUT_FILE.exists() and OUTPUT_FILE.stat().st_size < 50:
        OUTPUT_FILE.unlink()

    if OUTPUT_FILE.exists():
        print(f"Already exists: {OUTPUT_FILE}")
        return

    print("Querying GDC API for TCGA barcode-to-cancer-type mapping...")
    barcode_to_cancer: dict[str, str] = {}
    offset = 0

    while offset < MAX_CASES:
        print(f"  Fetching cases {offset}-{offset + BATCH_SIZE}...")
        try:
            data = fetch_gdc_cases(offset=offset, size=BATCH_SIZE)
        except Exception as e:
            print(f"  API error at offset {offset}: {e}")
            break

        hits = data.get("data", {}).get("hits", [])
        if not hits:
            break

        for hit in hits:
            submitter_id = hit.get("submitter_id", "")
            project = hit.get("project", {})
            project_id = project.get("project_id", "") if isinstance(project, dict) else ""

            if submitter_id and project_id:
                barcode_to_cancer[submitter_id] = project_id

        offset += BATCH_SIZE
        total = data.get("data", {}).get("pagination", {}).get("total", 0)
        print(f"  Got {len(hits)} hits (total: {total})")

        if offset >= total:
            break

        time.sleep(0.5)  # Rate-limit

    print(f"\nTotal mappings: {len(barcode_to_cancer)}")

    if not barcode_to_cancer:
        print("ERROR: No mappings fetched. Check network/API access.")
        return

    # Save as TSV
    with open(OUTPUT_FILE, "w") as f:
        f.write("barcode\tcancer_type\n")
        for barcode in sorted(barcode_to_cancer.keys()):
            f.write(f"{barcode}\t{barcode_to_cancer[barcode]}\n")

    print(f"Saved: {OUTPUT_FILE}")

    # Print summary
    cancer_counts: dict[str, int] = {}
    for cancer in barcode_to_cancer.values():
        cancer_counts[cancer] = cancer_counts.get(cancer, 0) + 1
    print("\nCancer type distribution:")
    for cancer, count in sorted(cancer_counts.items(), key=lambda x: -x[1]):
        print(f"  {cancer}: {count}")


if __name__ == "__main__":
    main()
