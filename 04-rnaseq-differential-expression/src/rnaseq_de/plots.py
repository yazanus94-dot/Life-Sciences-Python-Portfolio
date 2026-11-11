"""Figures."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


def volcano(res, path, alpha=0.05, lfc=1.0, label_genes=()):
    d = res.dropna(subset=["padj"]).copy()
    d["nl10"] = -np.log10(d["pvalue"].clip(lower=1e-300))
    sig = (d.padj < alpha) & (d.log2FC.abs() >= lfc)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(d.log2FC[~sig], d.nl10[~sig], s=4, c="lightgrey")
    ax.scatter(d.log2FC[sig], d.nl10[sig], s=6, c="crimson")
    for g in label_genes:
        if g in d.index:
            ax.annotate(g, (d.log2FC[g], d.nl10[g]), fontsize=8)
    ax.axvline(lfc, ls=":", c="k", lw=.6); ax.axvline(-lfc, ls=":", c="k", lw=.6)
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 p")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def ma_plot(res, path, alpha=0.05):
    d = res.dropna(subset=["padj"])
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    s = d.padj < alpha
    ax.scatter(np.log10(d.baseMean[~s] + 1), d.log2FC[~s], s=4, c="lightgrey")
    ax.scatter(np.log10(d.baseMean[s] + 1), d.log2FC[s], s=6, c="crimson")
    ax.axhline(0, c="k", lw=.6); ax.set_xlabel("log10(mean normalised count + 1)"); ax.set_ylabel("log2 fold change")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def pca_plot(counts, size_factors, coldata, color, shape, path):
    logn = np.log2(counts.div(size_factors, axis=1) + 1)
    top = logn.var(axis=1).sort_values(ascending=False).index[:1000]
    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(logn.loc[top].T.to_numpy())
    df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=logn.columns).join(coldata)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sns.scatterplot(data=df, x="PC1", y="PC2", hue=color, style=shape, s=90, ax=ax)
    ev = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({ev[0]:.0%})"); ax.set_ylabel(f"PC2 ({ev[1]:.0%})")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def top_gene_heatmap(counts, size_factors, res, path, n=30):
    genes = res.dropna(subset=["padj"]).index[:n]
    z = np.log2(counts.loc[genes].div(size_factors, axis=1) + 1)
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
    g = sns.clustermap(z, cmap="vlag", center=0, figsize=(6, 8), yticklabels=True)
    g.savefig(path, dpi=150); plt.close("all")
