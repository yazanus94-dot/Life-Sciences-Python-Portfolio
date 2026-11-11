"""Oral-dose compartment models: closed form, ODE, and multiple dosing.

Units follow the data: dose in mg/kg, time in h, V/F in L/kg, CL/F in L/h/kg.
F (bioavailability) is not identifiable from oral data alone, so V and CL are
apparent (V/F, CL/F).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


def conc_1cmt(t, dose, ka, ke, V):
    """One-compartment, first-order absorption and elimination (Bateman function)."""
    t = np.asarray(t, float)
    if abs(ka - ke) < 1e-8 * max(ka, ke):  # limit as ka -> ke
        c = dose * ka * t * np.exp(-ke * t) / V
    else:
        c = dose * ka / (V * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t))
    return np.where(t < 0, 0.0, c)


def conc_ode(t_eval, dose, params, model="1cmt"):
    """Solve the model ODEs numerically (tight tolerances).

    1cmt params: (ka, ke, V). 2cmt params: (ka, k10, k12, k21, V).
    States: gut amount, central amount, [peripheral amount].
    """
    t_eval = np.asarray(t_eval, float)
    if model == "1cmt":
        ka, ke, V = params
        def rhs(_, y):
            return [-ka * y[0], ka * y[0] - ke * y[1]]
        y0 = [dose, 0.0]
    elif model == "2cmt":
        ka, k10, k12, k21, V = params
        def rhs(_, y):
            return [-ka * y[0],
                    ka * y[0] - (k10 + k12) * y[1] + k21 * y[2],
                    k12 * y[1] - k21 * y[2]]
        y0 = [dose, 0.0, 0.0]
    else:
        raise ValueError(model)
    sol = solve_ivp(rhs, (0.0, float(t_eval.max())), y0, t_eval=t_eval,
                    method="LSODA", rtol=1e-10, atol=1e-12)
    return sol.y[1] / V


def conc_multi_dose(t, dose, tau, n_doses, ka, ke, V):
    """Superposition of identical oral doses every tau hours (linear PK)."""
    t = np.asarray(t, float)
    c = np.zeros_like(t)
    for k in range(n_doses):
        c += conc_1cmt(t - k * tau, dose, ka, ke, V)
    return c


def conc_steady_state(t, dose, tau, ka, ke, V):
    """Analytic steady-state profile within a dosing interval (0 <= t < tau)."""
    t = np.asarray(t, float)
    return dose * ka / (V * (ka - ke)) * (
        np.exp(-ke * t) / (1 - np.exp(-ke * tau)) - np.exp(-ka * t) / (1 - np.exp(-ka * tau)))
