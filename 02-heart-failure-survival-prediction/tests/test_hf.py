import numpy as np
import pytest

from hf_survival import classify, data, survival
from hf_survival.config import FEATURES
from hf_survival.synth import make_cohort


@pytest.fixture(scope="module")
def df():
    return make_cohort(299, seed=0)


def test_schema_validation(df, tmp_path):
    p = tmp_path / "c.csv"; df.to_csv(p, index=False)
    assert len(data.load(str(p))) == 299
    df.drop(columns=["age"]).to_csv(p, index=False)
    with pytest.raises(ValueError):
        data.load(str(p))


def test_km_logrank_detects_ef_effect(df):
    fits, p = survival.km_by_ef(df)
    assert p < 0.05 and len(fits) == 2


def test_cox_signs_match_generating_model(df):
    cph, _ = survival.fit_cox(df)
    s = cph.summary
    assert s.loc["ejection_fraction", "coef"] < 0
    assert s.loc["serum_creatinine", "coef"] > 0
    assert s.loc["age", "coef"] > 0


def test_cv_cindex_above_chance(df):
    assert survival.cv_cindex(df, n_repeats=2).mean() > 0.6


def test_classifier_auc_above_chance_without_time(df):
    assert classify.cv_auc(df, FEATURES, n_repeats=2).mean() > 0.6


def test_follow_up_time_leaks(df):
    d = classify.leakage_demo(df)
    assert d["leaky"] > d["honest"] + 0.05


def test_pipeline_scaler_is_inside_cv():
    m = classify.make_model()
    assert [n for n, _ in m.steps] == ["log", "scale", "clf"]
