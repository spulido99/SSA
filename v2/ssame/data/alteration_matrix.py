"""Alteration matrix construction from MAF files and v1 .m2 format."""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import overload

import numpy as np
import pandas as pd

from ssame.types import AlterationMatrix

logger = logging.getLogger(__name__)

try:
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.csv as pa_csv

    _HAS_PYARROW = True
except ImportError:
    _HAS_PYARROW = False

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
    use_pyarrow: bool | None = None,
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
        use_pyarrow: If True, use PyArrow CSV reader (faster). If False,
            use pandas. If None (default), auto-detect (PyArrow for
            uncompressed, pandas for .gz).

    Returns:
        Binary AlterationMatrix (gene x sample).
    """
    if variant_classes is None:
        variant_classes = _NON_SILENT_VARIANTS

    # Auto-detect backend: PyArrow can't read .gz natively
    if use_pyarrow is None:
        is_compressed = str(path).endswith(".gz")
        use_pyarrow = _HAS_PYARROW and not is_compressed

    if use_pyarrow and _HAS_PYARROW:
        return _parse_maf_pyarrow(path, min_samples, variant_classes, sample_barcodes)

    return _parse_maf_pandas(path, min_samples, variant_classes, sample_barcodes, chunksize)


def _parse_maf_pyarrow(
    path: Path,
    min_samples: int,
    variant_classes: frozenset[str],
    sample_barcodes: set[str] | None,
) -> AlterationMatrix:
    """Fast MAF parser using PyArrow."""
    # Skip comment lines by reading with autogenerated column names,
    # then filtering. PyArrow doesn't support comment= directly.
    # Instead, read only the needed columns.
    read_opts = pa_csv.ReadOptions(
        autogenerate_column_names=False,
    )
    parse_opts = pa_csv.ParseOptions(delimiter="\t")
    convert_opts = pa_csv.ConvertOptions(
        include_columns=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"],
    )

    table = pa_csv.read_csv(path, read_options=read_opts, parse_options=parse_opts, convert_options=convert_opts)

    # Filter: variant class in non-silent set
    vc_col = table.column("Variant_Classification")
    mask = pc.is_in(vc_col, value_set=pa.array(list(variant_classes)))
    table = table.filter(mask)

    # Filter: sample barcodes
    if sample_barcodes is not None:
        if len(sample_barcodes) == 0:
            return _build_matrix({}, min_samples)
        bc_col = table.column("Tumor_Sample_Barcode")
        bc_mask = pc.is_in(bc_col, value_set=pa.array(list(sample_barcodes)))
        table = table.filter(bc_mask)

    # Filter comment rows: Hugo_Symbol starting with '#'
    gene_col = table.column("Hugo_Symbol")
    comment_mask = pc.invert(pc.starts_with(gene_col, pattern="#"))
    table = table.filter(comment_mask)

    # Build gene → sample sets
    genes = table.column("Hugo_Symbol").to_pylist()
    samples = table.column("Tumor_Sample_Barcode").to_pylist()

    gene_samples: dict[str, set[str]] = defaultdict(set)
    for gene, sample in zip(genes, samples):
        gene_samples[gene].add(sample)

    return _build_matrix(gene_samples, min_samples)


def _parse_maf_pandas(
    path: Path,
    min_samples: int,
    variant_classes: frozenset[str],
    sample_barcodes: set[str] | None,
    chunksize: int,
) -> AlterationMatrix:
    """Original pandas chunked MAF parser."""
    _REQUIRED_COLS = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"]

    gene_samples: dict[str, set[str]] = defaultdict(set)

    for chunk in pd.read_csv(
        path,
        sep="\t",
        comment="#",
        usecols=_REQUIRED_COLS,
        chunksize=chunksize,
    ):
        chunk = chunk[chunk["Variant_Classification"].isin(variant_classes)]

        if sample_barcodes is not None:
            chunk = chunk[chunk["Tumor_Sample_Barcode"].isin(sample_barcodes)]

        for gene, sample in zip(chunk["Hugo_Symbol"], chunk["Tumor_Sample_Barcode"]):
            gene_samples[gene].add(sample)

    return _build_matrix(gene_samples, min_samples)


def _build_matrix(
    gene_samples: dict[str, set[str]],
    min_samples: int,
) -> AlterationMatrix:
    """Build AlterationMatrix from gene→sample dict."""
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


def load_tcga_cancer_map(path: Path) -> dict[str, str]:
    """Load TCGA barcode → cancer type mapping.

    Args:
        path: Path to TSV with columns 'barcode' and 'cancer_type'.

    Returns:
        Dict mapping TCGA barcode to cancer type (e.g. 'TCGA-BRCA').
    """
    if not path.exists():
        raise FileNotFoundError(f"TCGA cancer map not found: {path}")

    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["barcode"], df["cancer_type"]))
