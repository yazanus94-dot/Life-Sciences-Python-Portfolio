"""Run the validation experiments. Usage: python scripts/run_pipeline.py [--reps 20000]"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wfsim import analytic, plots, sim  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20000)
    ap.add_argument("--N", type=int, default=50)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    out = Path(a.out); (out / "figures").mkdir(parents=True, exist_ok=True)
    N, p0 = a.N, 1 / (2 * a.N)

    rows = []
    for i, s in enumerate([0.0, 0.01, 0.02, 0.05, 0.1]):
        fixed, t_abs = sim.fixation_experiment(a.reps, N, s, seed=42 + i)
        k = int(fixed.sum()); lo, hi = analytic.clopper_pearson(k, a.reps)
        th = analytic.kimura_fixation_prob(p0, N, s)
        rows.append(dict(s=s, fixations=k, reps=a.reps, p_hat=k / a.reps, lo=lo, hi=hi, kimura=th,
                         theory_in_ci=bool(lo <= th <= hi),
                         mean_t_fix=float(t_abs[fixed].mean()) if k else np.nan))
    tab = pd.DataFrame(rows); tab.to_csv(out / "fixation_probability.csv", index=False)
    print(f"Fixation probability of a new mutant, N={N}, {a.reps} replicates per s:")
    print(tab.round(4).to_string(index=False))
    print(f"\nNeutral theory: P(fix) = 1/(2N) = {p0:.4f}; mean time to fixation | fixed = "
          f"{analytic.neutral_time_to_fixation(p0, N):.1f} generations (simulated {tab.mean_t_fix[0]:.1f})")

    G = 100
    Nc = np.full(G, N); Nb = Nc.copy(); Nb[30:50] = max(N // 10, 2)
    Hc = (lambda t: 2 * t * (1 - t))(sim.trajectories(4000, Nc, 0.0, 0.5, seed=1)).mean(axis=1)
    Hb = (lambda t: 2 * t * (1 - t))(sim.trajectories(4000, Nb, 0.0, 0.5, seed=2)).mean(axis=1)
    Ec, Eb = analytic.expected_heterozygosity(0.5, Nc), analytic.expected_heterozygosity(0.5, Nb)
    print(f"Heterozygosity at generation {G}: constant N sim {Hc[-1]:.4f} vs theory {Ec[-1]:.4f}; "
          f"bottleneck sim {Hb[-1]:.4f} vs theory {Eb[-1]:.4f}")

    plots.fixation_vs_s(rows, N, out / "figures" / "fixation_vs_s.png")
    plots.heterozygosity(Hc, Hb, Nc, Nb, out / "figures" / "heterozygosity.png")
    plots.trajectory_panel(sim.trajectories(30, np.full(200, 100), 0.0, 0.1, seed=3),
                           sim.trajectories(30, np.full(200, 100), 0.05, 0.1, seed=3), 0.05,
                           out / "figures" / "trajectories.png")


if __name__ == "__main__":
    main()
