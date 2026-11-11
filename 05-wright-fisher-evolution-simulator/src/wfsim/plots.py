from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .analytic import expected_heterozygosity, kimura_fixation_prob


def trajectory_panel(traj_neutral, traj_sel, s, path):
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for a, tr, title in [(ax[0], traj_neutral, "Neutral (s = 0)"), (ax[1], traj_sel, f"Selection (s = {s})")]:
        a.plot(tr, lw=0.7, alpha=0.7); a.set_title(title); a.set_xlabel("Generation")
    ax[0].set_ylabel("Allele frequency")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def fixation_vs_s(rows, N, path):
    s = np.array([r["s"] for r in rows]); est = np.array([r["p_hat"] for r in rows])
    lo = np.array([r["lo"] for r in rows]); hi = np.array([r["hi"] for r in rows])
    ss = np.linspace(0, s.max() * 1.05, 200)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(ss, [kimura_fixation_prob(1 / (2 * N), N, x) for x in ss], "k-", label="Kimura")
    ax.plot(ss, 2 * ss + 0 * ss, "k:", lw=0.8, label="2s approximation")
    ax.errorbar(s, est, yerr=[est - lo, hi - est], fmt="o", capsize=3, label="simulation (95% CI)")
    ax.axhline(1 / (2 * N), c="grey", ls="--", lw=0.8, label="neutral 1/(2N)")
    ax.set_xlabel("selection coefficient s"); ax.set_ylabel("P(fixation of a new mutant)"); ax.legend(fontsize=8)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def heterozygosity(H_sim_const, H_sim_bn, N_const, N_bn, path):
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    g = np.arange(len(H_sim_const))
    ax.plot(g, H_sim_const, label="simulated, constant N"); ax.plot(g, expected_heterozygosity(0.5, N_const), "k--", lw=0.8)
    ax.plot(g, H_sim_bn, label="simulated, bottleneck"); ax.plot(g, expected_heterozygosity(0.5, N_bn), "k--", lw=0.8, label="theory")
    ax.set_xlabel("Generation"); ax.set_ylabel("Mean heterozygosity 2p(1-p)"); ax.legend(fontsize=8)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
