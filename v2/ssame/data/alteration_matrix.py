"""Alteration matrix construction from MAF files and v1 .m2 format."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import overload

import numpy as np
import pandas as pd

from ssame.types import AlterationMatrix

# Annotation pattern: gene_name(A) or gene_name(D) or locus(A/D)
_ANNOTATION_RE = re.compile(r"^(.+)\(([AD])\)$")

# MAF variant classifications considered non-silent
_NON_SILENT_VARIANTS = frozenset({
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Splice_Region",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Translation_Start_Site",
    "Nonstop_Mutation",
    "De_novo_Start_InFrame",
    "De_novo_Start_OutOfFrame",
})


@overload
def load_m2(
    path: Path, *, return_annotations: bool = False
) -> AlterationMatrix: ...


@overload
def load_m2(
    path: Path, *, return_annotations: bool = True
) -> tuple[AlterationMatrix, dict[str, dict[str, str]]]: ...


def load_m2(
    path: Path, *, return_annotations: bool = False
) -> AlterationMatrix | tuple[AlterationMatrix, dict[str, dict[str, str]]]:
    """Load a v1-format .m2 alteration matrix file.

    Each line is: sample_id<TAB>gene1<TAB>gene2<TAB>...
    Lines starting with # are comments.
    Genes may have annotations: GENE(A) = amplification, GENE(D) = deletion.

    Args:
        path: Path to the .m2 file.
        return_annotations: If True, also return a dict mapping
            (gene, sample) -> alteration_type.

    Returns:
        AlterationMatrix, optionally with annotations dict.
    """
    sample_genes: dict[str, set[str]] = {}
    # annotations: gene -> {sample -> type}
    annotations: dict[str, dict[str, str]] = defaultdict(dict)

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            sample_id = parts[0]
            genes_raw = parts[1:]
            clean_genes: set[str] = set()
            for raw_gene in genes_raw:
                match = _ANNOTATION_RE.match(raw_gene)
                if match:
                    gene_name = match.group(1)
                    alt_type = match.group(2)
                    annotations[gene_name][sample_id] = alt_type
                else:
                    gene_name = raw_gene
                clean_genes.add(gene_name)
            # Also add annotated genes to the gene set
            clean_genes.update(annotations.keys() & {
                m.group(1) for g in genes_raw if (m := _ANNOTATION_RE.match(g))
            })
            sample_genes[sample_id] = clean_genes

    # Build sorted gene and sample lists
    all_genes = sorted({g for genes in sample_genes.values() for g in genes})
    all_samples = sorted(sample_genes.keys())

    gene_index = {g: i for i, g in enumerate(all_genes)}
    data = np.zeros((len(all_genes), len(all_samples)), dtype=np.int8)

    for j, sample in enumerate(all_samples):
        for gene in sample_genes[sample]:
            if gene in gene_index:
                data[gene_index[gene], j] = 1

    matrix = AlterationMatrix(data=data, genes=all_genes, samples=all_samples)

    if return_annotations:
        return matrix, dict(annotations)
    return matrix


def export_m2(matrix: AlterationMatrix, path: Path) -> None:
    """Export an AlterationMatrix to v1-format .m2 file."""
    with open(path, "w") as f:
        for j, sample in enumerate(matrix.samples):
            altered_genes = [
                matrix.genes[i]
                for i in range(matrix.n_genes)
                if matrix.data[i, j] == 1
            ]
            if altered_genes:
                f.write(f"{sample}\t" + "\t".join(altered_genes) + "\n")


def parse_maf_to_matrix(
    path: Path,
    *,
    min_samples: int = 1,
    variant_classes: frozenset[str] | None = None,
    sample_barcodes: set[str] | None = None,
    chunksize: int = 500_000,
) -> AlterationMatrix:
    """Parse a MAF file into a binary alteration matrix.

    Reads in chunks to handle large files (e.g. MC3 pan-cancer MAF ~3GB).
    Supports .gz compressed files natively via pandas.

    Args:
        path: Path to MAF file (tab-separated, optionally .gz compressed).
        min_samples: Minimum number of samples a gene must be mutated in.
        variant_classes: Set of variant classifications to include.
            Defaults to non-silent mutations.
        sample_barcodes: Optional set of sample barcodes to include.
            If provided, only mutations in these samples are kept.
        chunksize: Number of rows per chunk for reading large files.

    Returns:
        Binary AlterationMatrix (gene x sample).
    """
    if variant_classes is None:
        variant_classes = _NON_SILENT_VARIANTS

    _REQUIRED_COLS = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"]

    gene_samples: dict[str, set[str]] = defaultdict(set)

    for chunk in pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=_REQUIRED_COLS,
        chunksize=chunksize,
    ):
        # Filter to non-silent variants
        chunk = chunk[chunk["Variant_Classification"].isin(variant_classes)]

        # Optionally filter by sample barcodes
        if sample_barcodes is not None:
            chunk = chunk[chunk["Tumor_Sample_Barcode"].isin(sample_barcodes)]

        # Accumulate gene → sample sets
        for gene, sample in zip(chunk["Hugo_Symbol"], chunk["Tumor_Sample_Barcode"]):
            gene_samples[gene].add(sample)

    # Filter by minimum sample count
    gene_samples = {g: s for g, s in gene_samples.items() if len(s) >= min_samples}

    all_genes = sorted(gene_samples.keys())
    all_samples = sorted({s for samples in gene_samples.values() for s in samples})

    gene_index = {g: i for i, g in enumerate(all_genes)}
    data = np.zeros((len(all_genes), len(all_samples)), dtype=np.int8)

    sample_index = {s: j for j, s in enumerate(all_samples)}
    for gene, samples in gene_samples.items():
        i = gene_index[gene]
        for sample in samples:
            data[i, sample_index[sample]] = 1

    return AlterationMatrix(data=data, genes=all_genes, samples=all_samples)
