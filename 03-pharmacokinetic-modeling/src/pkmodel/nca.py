"""Model-independent (non-compartmental) summaries, used to cross-check the model."""
from __future__ import annotations

import numpy as np


def nca(t, c, n_terminal: int = 4):
    """Cmax, Tmax, AUC0-tlast (linear trapezoid), terminal rate, AUC0-inf."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    auc_last = float(np.trapezoid(c, t))
    i = int(np.argmax(c))
    tt, cc = t[-n_terminal:], c[-n_terminal:]
    ok = cc > 0
    lam = -np.polyfit(tt[ok], np.log(cc[ok]), 1)[0]
    return dict(cmax=float(c[i]), tmax=float(t[i]), auc_last=auc_last, lambda_z=float(lam),
                auc_inf=auc_last + float(c[-1] / lam))
