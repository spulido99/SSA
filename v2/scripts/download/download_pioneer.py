"""Download PIONEER oncoPPI structural interaction predictions.

Source: Xiong et al., Nature Biotechnology 2024 (doi:10.1038/s41587-024-02428-4)
Data: https://pioneer.yulab.org/downloads

Downloads very_high confidence human predictions (166 MB).
"""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlretrieve

# PIONEER prediction files at different confidence levels
PIONEER_URLS = {
    "very_high": "https://pioneer.yulab.org/static/predictions/very_high/human.txt",
    "high": "https://pioneer.yulab.org/static/predictions/high/human.txt",
    "medium": "https://pioneer.yulab.org/static/predictions/medium/human.txt",
}

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "raw" / "networks" / "pioneer"


def download_pioneer(
    output_dir: Path | None = None,
    confidence_level: str = "very_high",
) -> Path:
    """Download PIONEER oncoPPI predictions.

    Args:
        output_dir: Target directory (default: data/raw/networks/pioneer/).
        confidence_level: One of 'very_high', 'high', 'medium'.

    Returns:
        Path to the downloaded prediction file.
    """
    out = output_dir or DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    if confidence_level not in PIONEER_URLS:
        msg = f"Unknown confidence level: {confidence_level}. Choose from: {list(PIONEER_URLS.keys())}"
        raise ValueError(msg)

    url = PIONEER_URLS[confidence_level]
    filename = f"pioneer_{confidence_level}_human.txt"
    target = out / filename

    if target.exists():
        print(f"Already exists: {target}")
        return target

    print(f"Downloading PIONEER {confidence_level} confidence predictions...")
    print(f"  URL: {url}")
    print(f"  Size: ~166 MB (very_high)")
    urlretrieve(url, target)
    print(f"  Saved: {target}")

    return target


if __name__ == "__main__":
    result = download_pioneer()
    print(f"\nDone. Data at: {result}")
