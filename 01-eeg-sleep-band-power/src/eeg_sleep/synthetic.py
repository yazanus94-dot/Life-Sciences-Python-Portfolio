"""Synthetic single-channel recordings with a known dominant frequency per stage.

Used by the offline tests and by `run_pipeline.py --synthetic` (demo only; the numbers
it produces are NOT findings about real sleep)."""
from __future__ import annotations

import mne
import numpy as np

SF = 100.0
INJECT = {"Sleep stage W": 10.0, "Sleep stage 1": 7.0, "Sleep stage 2": 14.0,
          "Sleep stage 3": 2.0, "Sleep stage R": 6.0}


def synthetic_raw(seed: int, n_per_stage: int = 8):
    rng = np.random.default_rng(seed)
    stages = list(INJECT) * n_per_stage
    rng.shuffle(stages)
    sig = []
    for s in stages:
        t = np.arange(0, 30, 1 / SF)
        sig.append(30e-6 * np.sin(2 * np.pi * INJECT[s] * t) + 5e-6 * rng.standard_normal(t.size))
    raw = mne.io.RawArray(np.concatenate(sig)[None, :], mne.create_info(["EEG Fpz-Cz"], SF, "eeg"), verbose="ERROR")
    raw.set_annotations(mne.Annotations(np.arange(len(stages)) * 30.0, 30.0, stages))
    return raw
