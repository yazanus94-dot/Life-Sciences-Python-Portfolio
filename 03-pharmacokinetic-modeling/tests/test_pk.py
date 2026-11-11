import numpy as np
import pytest

from pkmodel import fit, models, nca, population

TT = np.array([0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24.0])


def test_ode_matches_closed_form():
    p = (1.5, 0.1, 0.45)
    np.testing.assert_allclose(models.conc_ode(TT, 4.0, p), models.conc_1cmt(TT, 4.0, *p), rtol=1e-6)


def test_ka_equals_ke_limit_is_continuous():
    a = models.conc_1cmt(TT, 4.0, 0.1, 0.1, 0.45)
    b = models.conc_1cmt(TT, 4.0, 0.1 + 1e-6, 0.1, 0.45)
    np.testing.assert_allclose(a, b, rtol=1e-4)


def test_auc_inf_equals_dose_over_clearance():
    ka, ke, V, D = 1.5, 0.1, 0.45, 4.0
    t = np.linspace(0, 400, 400001)
    auc = np.trapezoid(models.conc_1cmt(t, D, ka, ke, V), t)
    assert auc == pytest.approx(D / (ke * V), rel=1e-4)


def test_two_compartment_reduces_to_one_when_no_exchange():
    p1 = (1.5, 0.1, 0.45)
    p2 = (1.5, 0.1, 1e-12, 1e-12, 0.45)
    np.testing.assert_allclose(models.conc_ode(TT, 4.0, p2, "2cmt"), models.conc_ode(TT, 4.0, p1), rtol=1e-6)


def test_parameter_recovery_with_noise():
    rng = np.random.default_rng(0)
    ka, ke, V = 1.5, 0.1, 0.45
    c = models.conc_1cmt(TT, 4.0, ka, ke, V) * np.exp(rng.normal(0, 0.05, TT.size))
    r = fit.fit_subject(TT, c, 4.0)
    ka_h, ke_h, V_h = r["params"]
    assert ka_h == pytest.approx(ka, rel=0.25) and ke_h == pytest.approx(ke, rel=0.1) and V_h == pytest.approx(V, rel=0.1)


def test_flip_flop_gives_identical_curve():
    """Swapping ka and ke (and rescaling V) reproduces the same curve exactly."""
    ka, ke, V, D = 1.5, 0.1, 0.45, 4.0
    V_swap = V * ke / ka  # D*ka/(V(ka-ke)) = D*ke/(V'(ke-ka)) * (-1)... equal when V' = V*ke/ka
    a = models.conc_1cmt(TT, D, ka, ke, V)
    b = models.conc_1cmt(TT, D, ke, ka, V_swap)
    np.testing.assert_allclose(a, b, rtol=1e-9)


def test_superposition_converges_to_analytic_steady_state():
    ka, ke, V, D, tau = 1.5, 0.1, 0.45, 4.0, 12.0
    t = np.linspace(0, tau, 50)
    multi = models.conc_multi_dose(t + 60 * tau, D, tau, 61, ka, ke, V)  # t after dose 61
    np.testing.assert_allclose(multi, models.conc_steady_state(t, D, tau, ka, ke, V), rtol=1e-6)


def test_nca_terminal_rate_recovers_ke():
    c = models.conc_1cmt(TT, 4.0, 1.5, 0.1, 0.45)
    assert nca.nca(TT, c)["lambda_z"] == pytest.approx(0.1, rel=0.02)


def test_aicc_penalises_extra_parameters():
    assert fit.aicc(1.0, 11, 5) > fit.aicc(1.0, 11, 3)


def test_real_data_loader_and_fit():
    df = population.load_theoph("data/Theoph.csv")
    assert df.Subject.nunique() == 12
    t = population.fit_all(df[df.Subject.isin([1, 2])])
    assert (t.ke > 0).all() and t.t_half.between(3, 20).all()
