"""Render docs/demo.gif: allele-frequency trajectories, neutral vs selected."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wfsim import sim  # noqa: E402

G, N, s = 200, 100, 0.03
neu = sim.trajectories(30, np.full(G, N), 0.0, 0.1, seed=11)
sel = sim.trajectories(30, np.full(G, N), s, 0.1, seed=11)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
lines = []
for a, title in zip(ax, ["Neutral (s = 0)", f"Selection (s = {s})"]):
    a.set_xlim(0, G); a.set_ylim(0, 1); a.set_title(title); a.set_xlabel("Generation")
    lines.append([a.plot([], [], lw=0.8, alpha=0.7)[0] for _ in range(30)])
ax[0].set_ylabel("Allele frequency")
frames = 80
def update(i):
    k = max(2, int((G + 1) * (i + 1) / frames))
    for ls, tr in zip(lines, (neu, sel)):
        for j, ln in enumerate(ls):
            ln.set_data(np.arange(k), tr[:k, j])
    return [ln for ls in lines for ln in ls]
FuncAnimation(fig, update, frames=frames, blit=True).save(
    Path(__file__).resolve().parents[1] / "docs" / "demo.gif", writer=PillowWriter(fps=15))
print("wrote docs/demo.gif")
