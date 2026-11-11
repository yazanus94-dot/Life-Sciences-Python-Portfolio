"""Analytic expectations (Kimura diffusion theory) and binomial confidence intervals."""
from __future__ import annotations

import numpy as np
from scipy.stats import beta, norm


def kimura_fixation_prob(p0, N, s):
    """P(fixation) = (1 - exp(-4 N s p0)) / (1 - exp(-4 N s)); neutral limit is p0.
    Matches the simulation model: 2N gene copies, allele fitness 1+s vs 1."""
    if abs(s) < 1e-12:
        return float(p0)
    return float(-np.expm1(-4 * N * s * p0) / -np.expm1(-4 * N * s))


def neutral_time_to_fixation(p0, N):
    """Mean generations to fixation, conditional on fixation, neutral:
    -4N (1-p) ln(1-p) / p  (-> 4N as p -> 0)."""
    return float(-4 * N * (1 - p0) * np.log(1 - p0) / p0)


def neutral_time_to_loss(p0, N):
    """Mean generations to loss, conditional on loss, neutral: -4N p ln(p) / (1-p)."""
    return float(-4 * N * p0 * np.log(p0) / (1 - p0))


def expected_heterozygosity(H0, N_schedule):
    """E[H_t] = H0 * prod(1 - 1/(2 N_i)); decay is governed by the harmonic mean of N."""
    f = 1 - 1 / (2 * np.asarray(N_schedule, float))
    return H0 * np.concatenate([[1.0], np.cumprod(f)])


def wald_interval(k, n, conf=0.95):
    p = k / n
    z = norm.ppf(0.5 + conf / 2)
    h = z * np.sqrt(p * (1 - p) / n)
    return max(p - h, 0.0), min(p + h, 1.0)


def wilson_interval(k, n, conf=0.95):
    z = norm.ppf(0.5 + conf / 2)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def clopper_pearson(k, n, conf=0.95):
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)
