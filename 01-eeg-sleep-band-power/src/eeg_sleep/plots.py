"""Figures."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from .config import BANDS, STAGE_ORDER
from .spectral import psd


def hypnogram_and_spectrogram(X, y, sfreq, path, max_freq=30):
    """Whole-night view for one subject: hypnogram above a time-frequency map."""
    f, p = psd(X, sfreq)
    keep = f <= max_freq
    logp = 10 * np.log10(p[:, keep] * 1e12 + 1e-12)  # uV^2/Hz in dB
    hours = np.arange(len(y)) * 30 / 3600
    order = {0: 4, 4: 3, 1: 2, 2: 1, 3: 0}  # AASM plotting order: W top, N3 bottom
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True,
                           gridspec_kw={"height_ratios": [1, 2]})
    ax[0].step(hours, [order[int(s)] for s in y], where="post", color="k", lw=1)
    ax[0].set_yticks(range(5)); ax[0].set_yticklabels(["N3", "N2", "N1", "REM", "W"])
    ax[0].set_ylabel("Stage")
    im = ax[1].pcolormesh(hours, f[keep], logp.T, shading="auto", cmap="viridis")
    ax[1].set_ylabel("Frequency (Hz)"); ax[1].set_xlabel("Time (h)")
    fig.colorbar(im, ax=ax, label="PSD (dB, uV^2/Hz)", pad=0.01)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def mean_psd_by_stage(X, y, sfreq, path, names):
    f, p = psd(X, sfreq)
    keep = (f >= 0.5) & (f <= 30)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for code, name in names.items():
        m = y == code
        if m.any():
            ax.semilogy(f[keep], np.median(p[m][:, keep], axis=0), label=f"{name} (n={m.sum()})")
    for lo, hi in BANDS.values():
        ax.axvline(lo, color="grey", lw=0.5, ls=":")
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Median PSD (V^2/Hz)"); ax.legend()
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def relative_power_heatmap(subj, path):
    cols = [c for c in subj.columns if c.endswith("_rel")]
    tab = subj.groupby("stage")[cols].median().loc[STAGE_ORDER]
    tab.columns = [c.replace("_rel", "") for c in tab.columns]
    fig, ax = plt.subplots(figsize=(6, 3.8))
    sns.heatmap(tab, annot=True, fmt=".2f", cmap="rocket_r", ax=ax, cbar_kws={"label": "relative power"})
    ax.set_ylabel("")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def subject_boxplots(subj, path):
    cols = [c for c in subj.columns if c.endswith("_rel")]
    long = subj.melt(id_vars=["subject", "stage"], value_vars=cols, var_name="band", value_name="rel")
    long["band"] = long["band"].str.replace("_rel", "")
    g = sns.catplot(data=long, x="stage", y="rel", col="band", kind="box",
                    order=STAGE_ORDER, col_wrap=3, height=2.6, sharey=False)
    g.savefig(path, dpi=150); plt.close("all")
