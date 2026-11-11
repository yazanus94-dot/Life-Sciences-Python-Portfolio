# Demo video script (4-5 min, screen recording with your voice)

Speak in your own words. Do not read this. You can also show docs/demo.gif on screen.

**0:00 Question (20 s).** How fast is theophylline absorbed and eliminated, and how variable is that between people?

**0:20 Data (30 s).** Theophylline dataset (Boeckmann, Sheiner, Beal 1994; R datasets "Theoph"): 12 subjects, one oral dose of about 4 mg/kg, 11 samples over 24 h.

**0:50 Models (80 s).** src/pkmodel/models.py: closed-form Bateman equation versus ODE solution (show the test that they agree), one- versus two-compartment, why V and CL are "apparent" (F not identifiable). fit.py: log-scale least squares (proportional error), log-parameters for positivity, multi-start, AICc.

**2:10 Run it live (40 s).** pytest, then python scripts/run_pipeline.py, then open results/figures/subject_fits.png.

**2:50 Results (60 s).** Report ONLY what your run produced: geometric-mean ka, ke, V/F, CL/F, half-life with bootstrap CIs; AICc preference; model AUC versus non-compartmental AUC. Mention any subject with a poorly identified ka (subject 9 in my run) and why.

**3:50 Evaluator problem (60 s).** Run problem_01_solution.py. Two fits, identical SSE, half-life 7 h versus 0.4 h. Explain flip-flop and how you would resolve it.

**4:50 Limitations (20 s).** Per-subject fits not a nonlinear mixed-effects model, single dose, no covariates, 11 points limit the 2-compartment model.
