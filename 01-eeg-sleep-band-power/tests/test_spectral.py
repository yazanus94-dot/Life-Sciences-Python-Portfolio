import numpy as np
import pytest

from eeg_sleep.config import BANDS
from eeg_sleep.preprocess import reject_artifacts
from eeg_sleep.spectral import band_powers

SF = 100.0
T = np.arange(0, 30, 1 / SF)


def test_pure_sine_lands_in_correct_band():
    rng = np.random.default_rng(0)
    for name, f0 in {"delta": 2.0, "theta": 6.0, "alpha": 10.0, "sigma": 14.0, "beta": 20.0}.items():
        x = np.sin(2 * np.pi * f0 * T) + 0.01 * rng.standard_normal(T.size)
        bp = band_powers(x[None, :], SF)
        assert bp[name]["rel"][0] > 0.95, name


def test_relative_power_sums_to_one():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((5, T.size))
    total = sum(band_powers(X, SF)[b]["rel"] for b in BANDS)
    np.testing.assert_allclose(total, 1.0, atol=1e-9)


def test_relative_power_is_amplitude_invariant():
    rng = np.random.default_rng(2)
    X = rng.standard_normal((3, T.size))
    a = band_powers(X, SF)["alpha"]["rel"]
    b = band_powers(10 * X, SF)["alpha"]["rel"]
    np.testing.assert_allclose(a, b, rtol=1e-9)


def test_absolute_power_matches_variance():
    # Parseval: integrating the PSD over all frequencies recovers the variance.
    from eeg_sleep.spectral import psd
    rng = np.random.default_rng(3)
    x = rng.standard_normal((1, 60 * int(SF)))
    f, p = psd(x, SF)
    assert np.isclose(p.sum() * (f[1] - f[0]), x.var(), rtol=0.1)


def test_artifact_rejection_threshold():
    X = np.zeros((3, 100)); X[1, 0] = 400e-6; X[1, 1] = -400e-6  # 800 uV p-p
    X[2, 0] = 100e-6
    y = np.array([0, 1, 2])
    Xc, yc, nrej = reject_artifacts(X, y, ptp_uv=500)
    assert nrej == 1 and list(yc) == [0, 2]
