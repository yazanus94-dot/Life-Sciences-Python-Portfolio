"""Filtering, epoching by hypnogram, and amplitude-based artifact rejection."""
from __future__ import annotations

import mne
import numpy as np

from .config import ANNOTATION_MAP, EPOCH_SEC


def bandpass(raw: mne.io.BaseRaw, l_freq: float = 0.3, h_freq: float = 35.0) -> mne.io.BaseRaw:
    """Zero-phase FIR band-pass (MNE default firwin, no phase distortion).

    Zero-phase matters: a causal filter shifts slow waves in time and distorts
    waveform shape. h_freq must stay below Nyquist (50 Hz at 100 Hz sampling).
    """
    return raw.copy().filter(l_freq, h_freq, phase="zero", fir_design="firwin", verbose="ERROR")


def epoch_by_stage(raw: mne.io.BaseRaw, epoch_sec: float = EPOCH_SEC):
    """Cut the single-channel recording into 30-s epochs labelled by stage.

    Returns (X, y): X has shape (n_epochs, n_samples) in volts; y holds stage codes.
    """
    events, _ = mne.events_from_annotations(
        raw, event_id=ANNOTATION_MAP, chunk_duration=epoch_sec, verbose="ERROR"
    )
    sfreq = raw.info["sfreq"]
    n = int(round(epoch_sec * sfreq))
    data = raw.get_data()[0]
    starts = events[:, 0] - raw.first_samp
    keep = (starts >= 0) & (starts + n <= data.shape[0])
    starts, labels = starts[keep], events[keep, 2]
    X = np.stack([data[s : s + n] for s in starts])
    return X, labels


def reject_artifacts(X: np.ndarray, y: np.ndarray, ptp_uv: float = 500.0):
    """Drop epochs whose peak-to-peak amplitude exceeds ptp_uv microvolts.

    500 uV removes gross artifacts (electrode pops, movement) while keeping
    genuine N3 slow waves, which commonly reach 150-300 uV peak-to-peak.
    Returns (X_clean, y_clean, n_rejected).
    """
    ptp = np.ptp(X, axis=1) * 1e6
    good = ptp <= ptp_uv
    return X[good], y[good], int((~good).sum())
