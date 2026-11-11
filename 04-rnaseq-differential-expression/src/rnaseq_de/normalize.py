"""Filtering and size-factor normalisation."""
from __future__ import annotations

import numpy as np
import pandas as pd


def filter_low_counts(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 2) -> pd.DataFrame:
    """Keep genes with at least min_count reads in at least min_samples samples."""
    return counts[(counts >= min_count).sum(axis=1) >= min_samples]


def median_of_ratios(counts: pd.DataFrame) -> pd.Series:
    """DESeq2-style size factors.

    Each sample is scaled by the median, over genes, of its ratio to the gene
    geometric mean across samples. Genes with any zero count are excluded from
    the reference, which makes the factor robust to a few dominant genes.
    """
    pos = counts[(counts > 0).all(axis=1)]
    if len(pos) < 10:
        raise ValueError("Too few genes with all-positive counts to estimate size factors.")
    logs = np.log(pos.to_numpy(float))
    ref = logs.mean(axis=1, keepdims=True)
    sf = np.exp(np.median(logs - ref, axis=0))
    return pd.Series(sf, index=counts.columns, name="size_factor")


def library_size_factors(counts: pd.DataFrame) -> pd.Series:
    """Total-count scaling (CPM-style), kept only to demonstrate its failure mode."""
    tot = counts.sum(axis=0).astype(float)
    return tot / np.exp(np.log(tot).mean())
