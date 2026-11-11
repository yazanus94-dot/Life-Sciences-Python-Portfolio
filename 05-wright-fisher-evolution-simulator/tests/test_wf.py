import numpy as np
import pytest

from wfsim import analytic, sim

N = 50
REPS = 40000


def se(p, n):
    return np.sqrt(p * (1 - p) / n)


def test_neutral_fixation_probability_is_1_over_2N():
    fixed, _ = sim.fixation_experiment(REPS, N, 0.0, seed=1)
    p0 = 1 / (2 * N)
    assert abs(fixed.mean() - p0) < 4 * se(p0, REPS)


def test_selected_fixation_probability_matches_kimura():
    s = 0.05
    fixed, _ = sim.fixation_experiment(REPS, N, s, seed=2)
    th = analytic.kimura_fixation_prob(1 / (2 * N), N, s)
    assert abs(fixed.mean() - th) < 4 * se(th, REPS)


def test_neutral_conditional_time_to_fixation():
    fixed, t = sim.fixation_experiment(REPS, N, 0.0, seed=3)
    th = analytic.neutral_time_to_fixation(1 / (2 * N), N)
    assert t[fixed].mean() == pytest.approx(th, rel=0.1)


def test_neutral_conditional_time_to_loss():
    fixed, t = sim.fixation_experiment(REPS, N, 0.0, seed=4)
    th = analytic.neutral_time_to_loss(1 / (2 * N), N)
    assert t[~fixed].mean() == pytest.approx(th, rel=0.1)


def _het(tr):
    return 2 * tr * (1 - tr)


def test_heterozygosity_decay_constant_N():
    G = 60
    h = _het(sim.trajectories(20000, np.full(G, N), 0.0, 0.5, seed=5))
    th = analytic.expected_heterozygosity(0.5, np.full(G, N))
    err = h.std(axis=1, ddof=1) / np.sqrt(h.shape[1])
    assert np.all(np.abs(h.mean(axis=1) - th) <= 4 * err + 1e-12)


def test_heterozygosity_decay_with_bottleneck():
    Ns = np.full(80, N); Ns[20:35] = 5
    h = _het(sim.trajectories(20000, Ns, 0.0, 0.5, seed=6))
    th = analytic.expected_heterozygosity(0.5, Ns)
    err = h.std(axis=1, ddof=1) / np.sqrt(h.shape[1])
    assert np.all(np.abs(h.mean(axis=1) - th) <= 4 * err + 1e-12)


def test_kimura_limits():
    assert analytic.kimura_fixation_prob(0.01, 100, 0.0) == 0.01
    assert analytic.kimura_fixation_prob(1 / 2000, 1000, 0.01) == pytest.approx(0.02, rel=0.02)
    assert analytic.kimura_fixation_prob(0.01, 100, -0.05) < 0.01  # deleterious


def test_boundaries_are_absorbing():
    tr = sim.trajectories(500, np.full(300, 20), 0.0, 0.2, seed=7)
    for j in range(tr.shape[1]):
        col = tr[:, j]
        hit = np.where((col == 0) | (col == 1))[0]
        if hit.size:
            assert np.all(col[hit[0]:] == col[hit[0]])


def test_seed_reproducibility():
    a, ta = sim.fixation_experiment(2000, 30, 0.02, seed=9)
    b, tb = sim.fixation_experiment(2000, 30, 0.02, seed=9)
    assert (a == b).all() and (ta == tb).all()


def test_interval_edge_cases():
    assert analytic.wald_interval(0, 100) == (0.0, 0.0)  # degenerate: the known Wald failure
    lo, hi = analytic.wilson_interval(0, 100)
    assert lo == pytest.approx(0, abs=1e-12) and hi > 0.02
    lo, hi = analytic.clopper_pearson(0, 100)
    assert lo == 0.0 and hi > 0.02
