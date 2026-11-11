from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .models import conc_1cmt, conc_multi_dose


def subject_fits(df, table, path):
    ids = sorted(df.Subject.unique())
    fig, axes = plt.subplots(3, 4, figsize=(12, 7.5), sharex=True)
    for ax, sid in zip(axes.ravel(), ids):
        g = df[df.Subject == sid]
        r = table[table.subject == sid].iloc[0]
        tt = np.linspace(0, g.Time.max(), 200)
        ax.plot(g.Time, g.conc, "ko", ms=4)
        ax.plot(tt, conc_1cmt(tt, r.dose, r.ka, r.ke, r.V_F), "C0-")
        ax.set_title(f"Subject {sid}  ka={r.ka:.2f} ke={r.ke:.3f}", fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Time (h)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Conc (mg/L)")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def accumulation(dose, tau, ka, ke, V, n_doses, path):
    t = np.linspace(0, tau * n_doses, 800)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, conc_multi_dose(t, dose, tau, n_doses, ka, ke, V))
    ax.set_xlabel("Time (h)"); ax.set_ylabel("Conc (mg/L)")
    ax.set_title(f"Repeated dosing every {tau:g} h (population geometric-mean parameters)")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
