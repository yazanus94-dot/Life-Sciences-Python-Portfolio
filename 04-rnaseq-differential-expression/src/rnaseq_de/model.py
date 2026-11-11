"""Per-gene negative-binomial GLM with Wald test and Benjamini-Hochberg FDR.

Design: log(mu) = log(size_factor) + X beta. X is built from coldata, e.g.
~ cell_line + dex for a paired design. The tested coefficient is the contrast
of interest (treated vs control).

Dispersion: per-gene Poisson-residual moment estimate (Cameron and Trivedi),
shrunk strongly (weight 0.1 on the gene estimate, log scale) toward a mean-dependent trend. Unshrunk gene-wise
estimates with ~8 samples were anti-conservative in my simulations (about 9
percent of null genes called at FDR 5 percent). P-values use a normal reference, as DESeq2 does (a t reference with
3 residual df was far too conservative here). This is a simplified stand-in for DESeq2; cross-check
borderline genes against PyDESeq2.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def design_matrix(coldata: pd.DataFrame, factors: list[str], contrast: str, ref: str) -> tuple[pd.DataFrame, str]:
    """Treatment-coded design. contrast is the factor of interest; ref its reference level.
    Returns (X, name of the coefficient to test)."""
    X = pd.get_dummies(coldata[factors].astype(str), drop_first=False, dtype=float)
    drop = [f"{contrast}_{ref}"]
    for f in factors:
        if f != contrast:
            first = sorted(coldata[f].astype(str).unique())[0]
            drop.append(f"{f}_{first}")
    X = X.drop(columns=drop)
    X.insert(0, "Intercept", 1.0)
    target = [c for c in X.columns if c.startswith(f"{contrast}_")]
    if len(target) != 1:
        raise ValueError("contrast must have exactly two levels")
    return X, target[0]


def _gene_dispersions(y_all, off, Xv):
    """Pass 1: per-gene moment dispersion from a Poisson fit."""
    n, p = Xv.shape
    out = np.full(len(y_all), np.nan)
    for i, y in enumerate(y_all):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mu = sm.GLM(y, Xv, family=sm.families.Poisson(), offset=off).fit().mu
            out[i] = (((y - mu) ** 2 - mu) / mu**2).sum() / max(n - p, 1)
        except Exception:
            pass
    return out


def _shrink_to_trend(alpha_gene, mean_norm, weight=0.1, floor=1e-4):
    """Shrink gene-wise dispersions toward a mean-dependent trend a0 + a1/mean
    (fit by least squares), averaging on the log scale. A simplified stand-in for
    DESeq2's empirical-Bayes shrinkage."""
    ok = np.isfinite(alpha_gene)
    a = np.clip(alpha_gene, floor, None)
    keep = ok & (a < np.nanpercentile(a[ok], 99))
    A = np.column_stack([np.ones(keep.sum()), 1.0 / mean_norm[keep]])
    coef = np.clip(np.linalg.lstsq(A, a[keep], rcond=None)[0], floor, None)
    trend = coef[0] + coef[1] / mean_norm
    return np.where(ok, np.exp(weight * np.log(a) + (1 - weight) * np.log(trend)), trend)


def fit_genes(counts: pd.DataFrame, size_factors: pd.Series, X: pd.DataFrame, target: str) -> pd.DataFrame:
    y_all = counts[X.index].to_numpy(float)
    sf = size_factors[X.index].to_numpy(float)
    off = np.log(sf)
    Xv = X.to_numpy(float)
    k = list(X.columns).index(target)
    n, p = Xv.shape
    mean_norm = (y_all / sf).mean(axis=1)
    alpha = _shrink_to_trend(_gene_dispersions(y_all, off, Xv), mean_norm)
    rows = []
    for g, y, a, m in zip(counts.index, y_all, alpha, mean_norm):
        beta = se = np.nan
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                nb = sm.GLM(y, Xv, family=sm.families.NegativeBinomial(alpha=float(a)), offset=off).fit()
                beta, se = float(nb.params[k]), float(nb.bse[k])
        except Exception:
            pass
        rows.append((g, float(m), beta / np.log(2), se / np.log(2), float(a)))
    res = pd.DataFrame(rows, columns=["gene", "baseMean", "log2FC", "lfcSE", "dispersion"]).set_index("gene")
    t = res["log2FC"] / res["lfcSE"]
    res["pvalue"] = 2 * norm.sf(np.abs(t))
    ok = res["pvalue"].notna()
    res["padj"] = np.nan
    res.loc[ok, "padj"] = multipletests(res.loc[ok, "pvalue"], method="fdr_bh")[1]
    return res.sort_values("padj")
