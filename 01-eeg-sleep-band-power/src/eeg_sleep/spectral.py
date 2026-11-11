"""Welch spectral estimation and band power."""
from __future__ import annotations

import numpy as np
from scipy.signal import welch

from .config import BANDS, TOTAL_RANGE


def psd(X: np.ndarray, sfreq: float, win_sec: float = 4.0):
    """Welch PSD per epoch: Hann window, 50 percent overlap, constant detrend.

    4-s windows give 0.25 Hz resolution, enough to separate delta from theta.
    """
    nperseg = int(win_sec * sfreq)
    return welch(
        X, fs=sfreq, window="hann", nperseg=nperseg, noverlap=nperseg // 2,
        detrend="constant", axis=-1,
    )


def band_powers(X: np.ndarray, sfreq: float, bands=BANDS, total=TOTAL_RANGE) -> dict:
    """Absolute and relative band power per epoch.

    Relative power = band power / power in total range, so values sum to 1 and
    are insensitive to subject-specific amplitude (skull thickness, electrodes).
    Returns dict band -> {"abs": array, "rel": array}.
    """
    f, p = psd(X, sfreq)
    df = f[1] - f[0]
    tot = p[:, (f >= total[0]) & (f < total[1])].sum(axis=1) * df
    out = {}
    for name, (lo, hi) in bands.items():
        m = (f >= lo) & (f < hi)
        a = p[:, m].sum(axis=1) * df
        out[name] = {"abs": a, "rel": a / tot}
    return out
