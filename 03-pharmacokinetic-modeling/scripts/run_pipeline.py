"""Fit PK models to the theophylline dataset.

Usage: python scripts/run_pipeline.py --csv data/Theoph.csv --out results
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pkmodel import plots, population  # noqa: E402

np.random.seed(42)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="data/Theoph.csv")
    ap.add_argument("--out", default="results")
    ap.add_argument("--tau", type=float, default=12.0, help="dosing interval for accumulation demo (h)")
    a = ap.parse_args()
    out = Path(a.out); (out / "figures").mkdir(parents=True, exist_ok=True)

    df = population.load_theoph(a.csv)
    table = population.fit_all(df)
    table.to_csv(out / "subject_fits.csv", index=False)
    print(table[["subject", "ka", "ke", "V_F", "CL_F", "t_half", "auc_inf_model", "auc_inf_nca"]].round(3).to_string(index=False))

    summ = []
    for p in ["ka", "ke", "V_F", "CL_F", "t_half"]:
        gm, lo, hi = population.bootstrap_geomean(table[p])
        summ.append(dict(parameter=p, geomean=gm, ci95_low=lo, ci95_high=hi))
    summ = pd.DataFrame(summ); summ.to_csv(out / "population_summary.csv", index=False)
    print("\nGeometric means (subject-level bootstrap 95% CI):\n", summ.round(3).to_string(index=False))

    d_aicc = table.aicc_2cmt - table.aicc_1cmt
    print(f"\nAICc: 1-compartment preferred (lower AICc) in {(d_aicc > 0).sum()} of {len(table)} subjects")
    ratio = table.auc_inf_model / table.auc_inf_nca
    print(f"AUC_inf model/NCA ratio: median {ratio.median():.3f}, range {ratio.min():.3f} to {ratio.max():.3f}")

    plots.subject_fits(df, table, out / "figures" / "subject_fits.png")
    g = {r.parameter: r.geomean for r in summ.itertuples()}
    plots.accumulation(float(table.dose.median()), a.tau, g["ka"], g["ke"], g["V_F"], 8,
                       out / "figures" / "accumulation.png")


if __name__ == "__main__":
    main()
