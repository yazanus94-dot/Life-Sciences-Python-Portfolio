"""Leakage-safe binary classification of DEATH_EVENT, plus a deliberate leakage demo.

Caveat: DEATH_EVENT within variable follow-up conflates risk with observation time,
so the Cox model is the primary analysis. The classifier is secondary."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

from .config import EVENT, FEATURES, SEED, TIME


def _logs(X: pd.DataFrame) -> pd.DataFrame:
    X = X.copy()
    for c in ("creatinine_phosphokinase", "serum_creatinine"):
        if c in X:
            X[c] = np.log(X[c])
    return X


def make_model(kind: str = "logreg") -> Pipeline:
    """Scaler lives INSIDE the pipeline so it is fit on training folds only."""
    clf = (LogisticRegression(max_iter=1000, C=1.0) if kind == "logreg"
           else RandomForestClassifier(n_estimators=300, min_samples_leaf=5, random_state=SEED))
    return Pipeline([("log", FunctionTransformer(_logs)), ("scale", StandardScaler()), ("clf", clf)])


def cv_auc(df: pd.DataFrame, features: list[str], kind: str = "logreg", n_repeats: int = 5) -> np.ndarray:
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=SEED)
    return cross_val_score(make_model(kind), df[features], df[EVENT], cv=cv, scoring="roc_auc")


def calibration(df: pd.DataFrame, features: list[str], kind: str = "logreg"):
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    p = cross_val_predict(make_model(kind), df[features], df[EVENT], cv=cv, method="predict_proba")[:, 1]
    frac, mean_pred = calibration_curve(df[EVENT], p, n_bins=5, strategy="quantile")
    return dict(frac_pos=frac, mean_pred=mean_pred, brier=brier_score_loss(df[EVENT], p),
                brier_baseline=brier_score_loss(df[EVENT], np.full(len(df), df[EVENT].mean())),
                auc=roc_auc_score(df[EVENT], p), p=p)


def leakage_demo(df: pd.DataFrame) -> dict:
    """Same model and CV; only difference is whether follow-up time is a feature."""
    return dict(honest=cv_auc(df, FEATURES).mean(), leaky=cv_auc(df, FEATURES + [TIME]).mean())
