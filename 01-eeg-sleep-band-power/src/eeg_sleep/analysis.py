"""Epoch table -> subject-level table -> paired statistics with FDR control."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon
from statsmodels.stats.multitest import multipletests

from .config import BANDS, STAGE_NAMES, STAGE_ORDER
from .spectral import band_powers


def epoch_table(X: np.ndarray, y: np.ndarray, sfreq: float, subject: int) -> pd.DataFrame:
    bp = band_powers(X, sfreq)
    df = pd.DataFrame({"subject": subject, "stage": [STAGE_NAMES[int(s)] for s in y]})
    for b in BANDS:
        df[f"{b}_rel"] = bp[b]["rel"]
        df[f"{b}_abs"] = bp[b]["abs"]
    return df


def subject_level(epochs: pd.DataFrame, value_suffix: str = "_rel") -> pd.DataFrame:
    """Median per subject x stage. The subject, not the epoch, is the unit of
    inference: epochs from one night are strongly autocorrelated."""
    cols = [f"{b}{value_suffix}" for b in BANDS]
    return epochs.groupby(["subject", "stage"])[cols].median().reset_index()


def stage_tests(subj: pd.DataFrame, ref: str = "W", suffix: str = "_rel") -> pd.DataFrame:
    """Per band: Friedman test across stages, plus paired Wilcoxon of each stage
    vs reference. Benjamini-Hochberg is applied across all stage-vs-ref tests
    (bands x stages) in one family. Only subjects with every stage are used."""
    rows = []
    for b in BANDS:
        wide = subj.pivot(index="subject", columns="stage", values=f"{b}{suffix}")
        wide = wide[STAGE_ORDER].dropna()
        n = len(wide)
        fr_p = friedmanchisquare(*[wide[s] for s in STAGE_ORDER]).pvalue if n >= 3 else np.nan
        for s in STAGE_ORDER:
            if s == ref:
                continue
            d = wide[s] - wide[ref]
            p = wilcoxon(d).pvalue if n >= 6 and (d != 0).any() else np.nan
            rows.append(dict(band=b, stage=s, ref=ref, n_subjects=n,
                             median_diff=d.median(), p_raw=p, friedman_p=fr_p))
    res = pd.DataFrame(rows)
    ok = res["p_raw"].notna()
    res.loc[ok, "p_fdr"] = multipletests(res.loc[ok, "p_raw"], method="fdr_bh")[1]
    return res
