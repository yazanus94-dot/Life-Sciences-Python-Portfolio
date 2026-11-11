"""Synthetic cohort with the same schema, used by tests and the evaluator problem.

Event times are exponential with a hazard driven by age, ejection fraction, serum
creatinine and sodium; follow-up is censored at a random administrative date, so the
observed time is informative about the event indicator (the leakage mechanism)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def make_cohort(n: int = 299, seed: int = 0) -> pd.DataFrame:
    r = np.random.default_rng(seed)
    d = pd.DataFrame({
        "age": r.normal(61, 12, n).clip(40, 95),
        "anaemia": r.integers(0, 2, n),
        "creatinine_phosphokinase": np.exp(r.normal(6.2, 1.0, n)).round(),
        "diabetes": r.integers(0, 2, n),
        "ejection_fraction": r.normal(38, 12, n).clip(14, 80).round(),
        "high_blood_pressure": r.integers(0, 2, n),
        "platelets": r.normal(263000, 97000, n).clip(25000, 850000),
        "serum_creatinine": np.exp(r.normal(0.25, 0.45, n)),
        "serum_sodium": r.normal(136.5, 4.4, n).clip(113, 148).round(),
        "sex": r.integers(0, 2, n),
        "smoking": r.integers(0, 2, n),
    })
    lin = (0.035 * (d.age - 60) - 0.05 * (d.ejection_fraction - 38)
           + 0.6 * (d.serum_creatinine - 1.3) - 0.04 * (d.serum_sodium - 136.5))
    t_event = r.exponential(1 / (0.0022 * np.exp(lin)))
    t_cens = r.uniform(60, 285, n)
    d["time"] = np.minimum(t_event, t_cens).round().clip(1)
    d["DEATH_EVENT"] = (t_event <= t_cens).astype(int)
    return d
