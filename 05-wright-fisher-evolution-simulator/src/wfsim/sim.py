"""Vectorised Wright-Fisher simulation.

Model: M = 2N gene copies. Each generation, selection maps frequency p to
p' = p(1+s) / (1+s p) (allele fitness 1+s vs 1), then drift draws the next count
from Binomial(M, p'). Frequencies 0 and 1 are absorbing.
"""
from __future__ import annotations

import numpy as np


def fixation_experiment(n_reps: int, N: int, s: float, p0: float | None = None,
                        seed: int = 42, max_gen: int = 10_000_000):
    """Run n_reps independent replicates until loss or fixation.

    Returns (fixed: bool array, t_abs: int array of absorption times).
    Default p0 = 1/(2N): a single new mutant copy.
    """
    rng = np.random.default_rng(seed)
    M = 2 * N
    k0 = 1 if p0 is None else int(round(p0 * M))
    k = np.full(n_reps, k0, dtype=np.int64)
    t_abs = np.zeros(n_reps, dtype=np.int64)
    active = np.arange(n_reps)
    t = 0
    while active.size and t < max_gen:
        t += 1
        p = k[active] / M
        p_sel = p * (1 + s) / (1 + s * p)
        k_new = rng.binomial(M, p_sel)
        k[active] = k_new
        done = (k_new == 0) | (k_new == M)
        t_abs[active[done]] = t
        active = active[~done]
    if active.size:
        raise RuntimeError(f"{active.size} replicates not absorbed by max_gen")
    return k == M, t_abs


def trajectories(n_reps: int, N_schedule, s: float = 0.0, p0: float = 0.5, seed: int = 42):
    """Allele-frequency trajectories for a (possibly time-varying) population size.

    N_schedule: array of diploid population sizes, one per generation.
    Returns array of shape (len(N_schedule) + 1, n_reps).
    """
    rng = np.random.default_rng(seed)
    Ns = np.asarray(N_schedule, dtype=np.int64)
    out = np.empty((len(Ns) + 1, n_reps))
    p = np.full(n_reps, p0)
    out[0] = p
    for g, N in enumerate(Ns):
        M = 2 * int(N)
        p_sel = p * (1 + s) / (1 + s * p)
        p = rng.binomial(M, p_sel) / M
        out[g + 1] = p
    return out
