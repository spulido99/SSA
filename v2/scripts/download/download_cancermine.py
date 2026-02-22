"""Download CancerMine text-mined cancer gene annotations.

Source: CancerMine (Zenodo, CC0 license)
Data: 856 drivers, 2,421 oncogenes, 2,037 TSGs across 426 cancer types
Paper: Nature Methods 2019 (doi:10.1038/s41592-019-0422-y)
"""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlretrieve

CANCERMINE_URL = "https://zenodo.org/records/16849846/files/cancermine_collated.tsv?download=1"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "raw" / "ground_truth"


def download_cancermine(output_dir: Path | None = None) -> Path:
    """Download CancerMine collated TSV.

    Returns:
        Path to cancermine_collated.tsv.
    """
    out = output_dir or DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    target = out / "cancermine_collated.tsv"
    if target.exists():
        print(f"Already exists: {target}")
        return target

    print("Downloading CancerMine from Zenodo...")
    print(f"  URL: {CANCERMINE_URL}")
    print("  Size: ~3.7 MB")
    urlretrieve(CANCERMINE_URL, target)
    print(f"  Saved: {target}")

    return target


if __name__ == "__main__":
    result = download_cancermine()
    print(f"\nDone. Data at: {result}")
