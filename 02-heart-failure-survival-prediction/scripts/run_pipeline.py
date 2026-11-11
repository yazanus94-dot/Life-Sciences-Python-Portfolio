"""Run the analysis. Usage: python scripts/run_pipeline.py --csv data/heart_failure_clinical_records_dataset.csv

If --csv is omitted the script tries to download from the UCI repository; if that fails,
download the CSV manually from the dataset page and pass its path."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hf_survival import classify, data, plots, survival  # noqa: E402
from hf_survival.config import FEATURES  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=data.URLS[0])
    ap.add_argument("--out", default="results")
    ap.add_argument("--synthetic", action="store_true", help="demo on a synthetic cohort (no download)")
    a = ap.parse_args()
    out = Path(a.out); (out / "figures").mkdir(parents=True, exist_ok=True)
    np.random.seed(42)

    if a.synthetic:
        from hf_survival.synth import make_cohort
        df = make_cohort(299, seed=0)
    else:
        df = data.load(a.csv)
    print(f"{len(df)} patients, {df.DEATH_EVENT.sum()} events ({df.DEATH_EVENT.mean():.1%})")

    fits, p = survival.km_by_ef(df)
    plots.km_plot(fits, p, out / "figures" / "km_ef.png")
    cph, ph = survival.fit_cox(df)
    cph.summary.to_csv(out / "cox_summary.csv"); ph.to_csv(out / "cox_ph_test.csv")
    plots.forest_plot(cph, out / "figures" / "cox_forest.png")
    print(cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3))
    print("PH test (rank):\n", ph.round(3))
    ci = survival.cv_cindex(df)
    print(f"Cox CV C-index: {ci.mean():.3f} +/- {ci.std():.3f}")

    rows = []
    for kind in ("logreg", "forest"):
        s = classify.cv_auc(df, FEATURES, kind)
        rows.append(dict(model=kind, auc_mean=s.mean(), auc_sd=s.std()))
    pd.DataFrame(rows).to_csv(out / "classifier_cv_auc.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    cal = classify.calibration(df, FEATURES)
    plots.calibration_plot(cal, out / "figures" / "calibration.png")
    print(f"Calibration: Brier {cal["brier"]:.3f} (baseline {cal["brier_baseline"]:.3f})")
    leak = classify.leakage_demo(df)
    pd.DataFrame([leak]).to_csv(out / "leakage_demo.csv", index=False)
    print(f"Leakage demo AUC: honest {leak["honest"]:.3f} vs with follow-up time {leak["leaky"]:.3f}")


if __name__ == "__main__":
    main()
