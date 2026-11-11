"""Kaplan-Meier, log-rank, Cox regression, PH check, cross-validated C-index."""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.utils import concordance_index
from sklearn.model_selection import RepeatedStratifiedKFold

from .config import EF_CUTOFF, EVENT, FEATURES, SEED, TIME

# Skewed positive covariates are log-transformed before the Cox model.
LOG_COLS = ["creatinine_phosphokinase", "serum_creatinine"]


def prep_cox(df: pd.DataFrame) -> pd.DataFrame:
    x = df[FEATURES + [TIME, EVENT]].copy()
    for c in LOG_COLS:
        x[c] = np.log(x[c])
    x["platelets"] = x["platelets"] / 1e5  # per 100k, readable hazard ratios
    return x


def km_by_ef(df: pd.DataFrame):
    low = df.ejection_fraction < EF_CUTOFF
    lr = logrank_test(df[TIME][low], df[TIME][~low], df[EVENT][low], df[EVENT][~low])
    fits = {}
    for name, m in {f"EF < {EF_CUTOFF}": low, f"EF >= {EF_CUTOFF}": ~low}.items():
        fits[name] = KaplanMeierFitter().fit(df[TIME][m], df[EVENT][m], label=f"{name} (n={m.sum()})")
    return fits, lr.p_value


def fit_cox(df: pd.DataFrame, penalizer: float = 0.01):
    x = prep_cox(df)
    mu, sd = x[FEATURES].mean(), x[FEATURES].std()
    x[FEATURES] = (x[FEATURES] - mu) / sd  # HR per 1 SD, comparable across covariates
    cph = CoxPHFitter(penalizer=penalizer).fit(x, TIME, EVENT)
    ph = proportional_hazard_test(cph, x, time_transform="rank")
    return cph, ph.summary


def cv_cindex(df: pd.DataFrame, n_splits: int = 5, n_repeats: int = 5, penalizer: float = 0.01) -> np.ndarray:
    """Out-of-fold C-index. Scaling uses training-fold statistics only."""
    x = prep_cox(df)
    scores = []
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)
    for tr, te in cv.split(x, x[EVENT]):
        a, b = x.iloc[tr].copy(), x.iloc[te].copy()
        mu, sd = a[FEATURES].mean(), a[FEATURES].std()
        a[FEATURES], b[FEATURES] = (a[FEATURES] - mu) / sd, (b[FEATURES] - mu) / sd
        m = CoxPHFitter(penalizer=penalizer).fit(a, TIME, EVENT)
        risk = m.predict_partial_hazard(b)
        scores.append(concordance_index(b[TIME], -risk, b[EVENT]))
    return np.array(scores)
