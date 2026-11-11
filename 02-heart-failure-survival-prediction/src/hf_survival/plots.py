from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def km_plot(fits, p, path):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for f in fits.values():
        f.plot_survival_function(ax=ax, ci_show=True)
    ax.set_xlabel("Follow-up (days)"); ax.set_ylabel("Survival probability")
    ax.set_title(f"Kaplan-Meier by ejection fraction (log-rank p = {p:.2g})")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def forest_plot(cph, path):
    s = cph.summary.sort_values("exp(coef)")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    y = np.arange(len(s))
    ax.errorbar(s["exp(coef)"], y, xerr=[s["exp(coef)"] - s["exp(coef) lower 95%"],
                s["exp(coef) upper 95%"] - s["exp(coef)"]], fmt="o", color="k", capsize=3)
    ax.axvline(1, ls=":", c="grey"); ax.set_yticks(y); ax.set_yticklabels(s.index)
    ax.set_xlabel("Hazard ratio per 1 SD (95% CI)")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def calibration_plot(cal, path):
    fig, ax = plt.subplots(figsize=(4.8, 4.5))
    ax.plot([0, 1], [0, 1], ":", c="grey"); ax.plot(cal["mean_pred"], cal["frac_pos"], "o-")
    ax.set_xlabel("Mean predicted risk"); ax.set_ylabel("Observed fraction")
    ax.set_title(f"Calibration (Brier {cal["brier"]:.3f} vs baseline {cal["brier_baseline"]:.3f})", fontsize=9)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
