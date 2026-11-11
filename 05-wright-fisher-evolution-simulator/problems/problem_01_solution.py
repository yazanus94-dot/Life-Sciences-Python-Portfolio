"""Reference solution for problems/problem_01_fixation_estimate.md.

N = 1000 diploids, s = 0.01, a single new beneficial mutant. Compare:
  (a) a 100-replicate estimate with a Wald interval (the common mistake),
  (b) the same with Wilson / Clopper-Pearson intervals,
  (c) an adequately powered run (200000 replicates) against Kimura.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wfsim import analytic, sim  # noqa: E402

N, s = 1000, 0.01
p0 = 1 / (2 * N)
u = analytic.kimura_fixation_prob(p0, N, s)
print(f"Kimura u = {u:.5f}   (2s approximation = {2 * s:.3f}, neutral 1/(2N) = {p0:.4f})\n")

# (a)+(b): coverage of 95% intervals when n = 100 replicates are used
rng = np.random.default_rng(42)
n, trials = 100, 20000
ks = rng.binomial(n, u, trials)
for name, f in [("Wald", analytic.wald_interval), ("Wilson", analytic.wilson_interval),
                ("Clopper-Pearson", analytic.clopper_pearson)]:
    cover = np.mean([f(k, n)[0] <= u <= f(k, n)[1] for k in ks])
    print(f"{name:16s} coverage with n=100: {cover:.3f}")
print(f"P(zero fixations in 100 replicates) = {(1 - u) ** n:.3f}  -> Wald interval collapses to [0, 0]\n")

# (c): a properly powered simulation
reps = 200_000
fixed, _ = sim.fixation_experiment(reps, N, s, seed=42)
k = int(fixed.sum()); lo, hi = analytic.clopper_pearson(k, reps)
print(f"Simulation with {reps} replicates: {k} fixations, estimate {k / reps:.5f}, 95% CI [{lo:.5f}, {hi:.5f}]")
print(f"Kimura {u:.5f} inside CI: {lo <= u <= hi}")
