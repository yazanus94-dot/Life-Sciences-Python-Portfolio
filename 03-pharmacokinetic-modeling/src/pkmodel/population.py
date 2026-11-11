"""Per-subject fits, summaries, and subject-level bootstrap."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .fit import fit_subject
from .nca import nca


def load_theoph(path: str) -> pd.DataFrame:
    d = pd.read_csv(path)
    need = {"Subject", "Dose", "Time", "conc"}
    if not need <= set(d.columns):
        raise ValueError(f"missing columns {need - set(d.columns)}")
    return d


def fit_all(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sid, g in df.groupby("Subject"):
        t, c, dose = g.Time.to_numpy(), g.conc.to_numpy(), float(g.Dose.iloc[0])
        f1 = fit_subject(t, c, dose, "1cmt")
        f2 = fit_subject(t, c, dose, "2cmt")
        ka, ke, V = f1["params"]
        nc = nca(t, c)
        rows.append(dict(subject=sid, dose=dose, ka=ka, ke=ke, V_F=V, CL_F=ke * V,
                         t_half=np.log(2) / ke, auc_inf_model=dose / (ke * V), auc_inf_nca=nc["auc_inf"],
                         cmax_obs=nc["cmax"], tmax_obs=nc["tmax"],
                         sse_1cmt=f1["sse"], aicc_1cmt=f1["aicc"], aicc_2cmt=f2["aicc"], n_obs=f1["n"]))
    return pd.DataFrame(rows)


def bootstrap_geomean(values, n_boot: int = 4000, seed: int = 42, ci: float = 95):
    """Percentile bootstrap CI for the geometric mean, resampling SUBJECTS."""
    v = np.log(np.asarray(values, float))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), (n_boot, len(v)))
    boots = np.exp(v[idx].mean(axis=1))
    lo, hi = np.percentile(boots, [(100 - ci) / 2, 100 - (100 - ci) / 2])
    return float(np.exp(v.mean())), float(lo), float(hi)
