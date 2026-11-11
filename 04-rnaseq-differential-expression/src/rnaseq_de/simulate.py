"""Synthetic paired RNA-seq counts with known truth (used by tests and the problem)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate(n_genes=800, n_de=80, lfc=2.0, n_pairs=4, seed=0, disp=0.1, dominant=None):
    """4 subjects x (control, treated). Returns counts, coldata, set of true DE genes.

    dominant: optional (n_genes_dominant, fold) -> a few very highly expressed genes are
    up-regulated in treated samples, which distorts total-count normalisation.
    """
    rng = np.random.default_rng(seed)
    base = np.exp(rng.normal(5.0, 1.5, n_genes))
    de = rng.choice(n_genes, n_de, replace=False)
    sign = rng.choice([-1, 1], n_de)
    cols, trt, subj = [], [], []
    depth = rng.uniform(0.7, 1.4, 2 * n_pairs)
    X = np.zeros((n_genes, 2 * n_pairs))
    j = 0
    for s in range(n_pairs):
        subj_eff = np.exp(rng.normal(0, 0.3, n_genes))
        for t in (0, 1):
            mu = base * subj_eff * depth[j]
            if t:
                mu = mu.copy(); mu[de] *= 2.0 ** (sign * lfc)
                if dominant:
                    k, fold = dominant
                    top = np.argsort(base)[-k:]
                    mu[top] *= fold
            X[:, j] = rng.negative_binomial(1 / disp, 1 / (1 + disp * mu))
            label = "trt" if t else "ctl"
            cols.append(f"S{s}_{label}"); trt.append(label); subj.append(f"P{s}")
            j += 1
    genes = [f"g{i}" for i in range(n_genes)]
    counts = pd.DataFrame(X.astype(int), index=genes, columns=cols)
    coldata = pd.DataFrame({"subject": subj, "treat": trt}, index=cols)
    truth = {genes[i] for i in de}
    return counts, coldata, truth
