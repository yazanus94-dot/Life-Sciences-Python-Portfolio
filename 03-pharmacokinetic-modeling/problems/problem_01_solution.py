"""Reference solution for problems/problem_01_flip_flop.md.

Simulate oral data (ka=1.5/h, ke=0.1/h, V/F=0.45 L/kg), then fit the same data from two
starting points: one with ka > ke, one with ka < ke. Both fits are equally good; only
one gives the right half-life.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pkmodel import fit, models  # noqa: E402

rng = np.random.default_rng(42)
t = np.array([0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24.0])
true = (1.5, 0.1, 0.45)
c = models.conc_1cmt(t, 4.0, *true) * np.exp(rng.normal(0, 0.05, t.size))

print(f"truth: ka={true[0]}, ke={true[1]}, V/F={true[2]}, half-life={np.log(2) / true[1]:.2f} h\n")
for label, start in [("start ka > ke", (1.5, 0.1, 0.45)), ("start ka < ke", (0.1, 1.5, 0.045))]:
    r = fit.fit_subject(t, c, 4.0, starts=[start])
    ka, ke, V = r["params"]
    print(f"{label}: ka={ka:.3f} ke={ke:.3f} V/F={V:.3f}  SSE={r['sse']:.5f}  'half-life' ln2/ke={np.log(2) / ke:.2f} h")
