"""Download all open-access ground truth databases.

Runs all three download scripts:
- IntOGen (~600 drivers)
- NCG (~3,347 drivers)
- CancerMine (~5,000 annotated genes)

Usage:
    python scripts/download/download_ground_truth.py
"""

from __future__ import annotations

from download_cancermine import download_cancermine
from download_intogen import download_intogen
from download_ncg import download_ncg


def main() -> None:
    print("=" * 60)
    print("Downloading open-access ground truth databases")
    print("=" * 60)

    print("\n[1/3] IntOGen")
    intogen = download_intogen()

    print("\n[2/3] NCG")
    ncg = download_ncg()

    print("\n[3/3] CancerMine")
    cancermine = download_cancermine()

    print("\n" + "=" * 60)
    print("Done! Files:")
    print(f"  IntOGen:    {intogen}")
    print(f"  NCG:        {ncg}")
    print(f"  CancerMine: {cancermine}")
    print("=" * 60)


if __name__ == "__main__":
    main()
