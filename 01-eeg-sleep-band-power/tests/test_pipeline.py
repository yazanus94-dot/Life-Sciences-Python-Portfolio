"""End-to-end test on a synthetic recording with a known ground truth,
so it runs offline without downloading PhysioNet data."""
import mne
import numpy as np
import pandas as pd

from eeg_sleep import analysis
from eeg_sleep.preprocess import bandpass, epoch_by_stage

from eeg_sleep.synthetic import INJECT, SF, synthetic_raw  # noqa: E402,F401

EXPECTED_BAND = {"W": "alpha", "N1": "theta", "N2": "sigma", "N3": "delta", "REM": "theta"}


def _subject_table():
    tables = []
    for sid in range(8):
        raw = bandpass(synthetic_raw(sid))
        X, y = epoch_by_stage(raw)
        tables.append(analysis.epoch_table(X, y, SF, sid))
    return analysis.subject_level(pd.concat(tables))


def test_recovers_injected_band_per_stage():
    subj = _subject_table()
    med = subj.groupby("stage")[[c for c in subj.columns if c.endswith("_rel")]].median()
    for stage, band in EXPECTED_BAND.items():
        assert med.loc[stage].idxmax() == f"{band}_rel", stage


def test_fdr_adjusted_p_never_smaller_than_raw():
    stats = analysis.stage_tests(_subject_table())
    ok = stats.dropna(subset=["p_raw"])
    assert len(ok) > 0
    assert (ok.p_fdr >= ok.p_raw - 1e-12).all()
