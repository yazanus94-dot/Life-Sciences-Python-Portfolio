"""Nonlinear least-squares fitting on log-transformed concentrations.

Fitting log(c) assumes proportional (multiplicative) error, which is typical for
assay data, and keeps all parameters positive via a log parameterisation.
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import least_squares

from .models import conc_1cmt, conc_ode

K = {"1cmt": 3, "2cmt": 5}


def aicc(sse: float, n: int, k: int) -> float:
    """AICc for least squares on log scale; k counts the error variance too."""
    kk = k + 1
    if n - kk - 1 <= 0:
        return np.inf
    return n * np.log(sse / n) + 2 * kk + 2 * kk * (kk + 1) / (n - kk - 1)


def _predict(theta, t, dose, model):
    p = np.exp(theta)
    return conc_1cmt(t, dose, *p) if model == "1cmt" else conc_ode(t, dose, p, "2cmt")


def fit_subject(t, c, dose, model="1cmt", starts=None):
    """Multi-start fit. Returns dict(params, sse, aicc, n, success).

    Observations with c <= 0 are dropped (log undefined)."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    m = (c > 0) & (t > 0)
    t, c = t[m], c[m]
    if starts is None:
        if model == "1cmt":
            starts = [(ka, ke, v) for ka, ke, v in itertools.product([0.5, 1.5, 4.0], [0.05, 0.1, 0.3], [0.3, 0.6])]
        else:
            starts = [(ka, k10, 0.3, 0.2, v) for ka, k10, v in itertools.product([0.5, 2.0], [0.1, 0.3], [0.3, 0.6])]
    best = None
    for s in starts:
        try:
            r = least_squares(lambda th: np.log(c) - np.log(np.maximum(_predict(th, t, dose, model), 1e-12)),
                              np.log(s), method="trf", x_scale=1.0, max_nfev=400)
        except Exception:
            continue
        sse = float(np.sum(r.fun ** 2))
        if best is None or sse < best["sse"]:
            best = dict(params=np.exp(r.x), sse=sse, success=bool(r.success))
    if best is None:
        raise RuntimeError("all starts failed")
    best["n"] = len(c)
    best["aicc"] = aicc(best["sse"], len(c), K[model])
    return best
