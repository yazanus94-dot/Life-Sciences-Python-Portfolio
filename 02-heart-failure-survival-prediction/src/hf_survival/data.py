"""Loading and validating the UCI Heart Failure Clinical Records CSV."""
from __future__ import annotations

import pandas as pd

from .config import EVENT, REQUIRED, TIME

URLS = [
    "https://archive.ics.uci.edu/ml/machine-learning-databases/00519/heart_failure_clinical_records_dataset.csv",
]


def load(path_or_url: str) -> pd.DataFrame:
    df = pd.read_csv(path_or_url)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if df[REQUIRED].isna().any().any():
        raise ValueError("Unexpected missing values")
    if not set(df[EVENT].unique()) <= {0, 1}:
        raise ValueError("DEATH_EVENT must be 0/1")
    if (df[TIME] <= 0).any():
        raise ValueError("time must be positive")
    return df
