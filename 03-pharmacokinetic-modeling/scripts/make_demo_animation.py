"""Render docs/demo.gif: drug accumulating to steady state under repeated dosing,
with the analytic steady-state profile overlaid. Embed in the README."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pkmodel.models import conc_multi_dose, conc_steady_state  # noqa: E402

ka, ke, V, D, tau, N = 1.36, 0.088, 0.454, 4.0, 12.0, 8  # approx. population geometric means
t = np.linspace(0, tau * N, 800)
c = conc_multi_dose(t, D, tau, N, ka, ke, V)
tt = np.linspace(0, tau, 100)
ss = conc_steady_state(tt, D, tau, ka, ke, V)

fig, ax = plt.subplots(figsize=(7, 3.8))
ax.set_xlim(0, tau * N); ax.set_ylim(0, c.max() * 1.1)
ax.set_xlabel("Time (h)"); ax.set_ylabel("Concentration (mg/L)")
ax.set_title("Accumulation to steady state (dose every 12 h)")
(line,) = ax.plot([], [], lw=2)
ax.plot((N - 1) * tau + tt, ss, "r--", lw=1, label="analytic steady state")
ax.legend(loc="lower right")
frames = 60
def update(i):
    k = max(2, int(len(t) * (i + 1) / frames))
    line.set_data(t[:k], c[:k]); return (line,)
FuncAnimation(fig, update, frames=frames, blit=True).save(
    Path(__file__).resolve().parents[1] / "docs" / "demo.gif", writer=PillowWriter(fps=15))
print("wrote docs/demo.gif")
