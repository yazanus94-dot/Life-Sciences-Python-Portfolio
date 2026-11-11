"""Reference solution for problems/problem_01_leakage.md.

Same model, same cross-validation; the only change is whether follow-up time is a
feature. Runs on a synthetic cohort by default (schema identical to the UCI file),
or on the real CSV if a path is given:  python problems/problem_01_solution.py data/heart_failure_clinical_records_dataset.csv
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hf_survival import classify, data  # noqa: E402
from hf_survival.config import FEATURES, TIME  # noqa: E402
from hf_survival.synth import make_cohort  # noqa: E402

if len(sys.argv) > 1:
    cohorts = {"real data": data.load(sys.argv[1])}
else:
    cohorts = {f"synthetic seed {s}": make_cohort(299, s) for s in range(5)}

for name, df in cohorts.items():
    honest = classify.cv_auc(df, FEATURES).mean()
    leaky = classify.cv_auc(df, FEATURES + [TIME]).mean()
    print(f"{name:20s} CV AUC without time: {honest:.3f}   with time: {leaky:.3f}   inflation: {leaky - honest:+.3f}")
