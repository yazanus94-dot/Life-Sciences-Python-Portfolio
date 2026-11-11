# Demo video script (4-5 min, screen recording with your voice)

Speak in your own words. Do not read this.

**0:00 Question (20 s).** Which routinely recorded clinical measurements are associated with death in heart-failure patients, and can I predict risk without cheating?

**0:20 Data (30 s).** UCI Heart Failure Clinical Records: 299 patients, 12 clinical features, follow-up time and death event. Mention it is small, from one hospital cohort, and what that means for generalisation.

**0:50 Survival analysis (70 s).** Kaplan-Meier by ejection fraction below or above 40 (pre-specified clinical cutoff, not tuned to the data), log-rank test, then multivariable Cox model. Explain hazard ratios per 1 SD and the proportional-hazards check.

**2:00 Prediction (60 s).** Logistic regression and random forest with repeated stratified 5-fold CV, scaler inside a Pipeline, calibration curve, Brier score against the prevalence baseline. State plainly that the Cox model is primary because follow-up differs by patient.

**3:00 Run it live (30 s).** pytest, then run_pipeline.py. Show results/ files.

**3:30 Results (50 s).** Report ONLY what your run produced: which covariates are significant, C-index, AUC, Brier. Compare with the published work (Chicco and Jurman 2020 highlight serum creatinine and ejection fraction) and say where you agree or differ.

**4:20 Evaluator problem (50 s).** Run problem_01_solution.py on the real CSV. Explain why follow-up time is leakage and show the AUC with and without it.

**5:10 Limitations (20 s).** Small n, single cohort, no external validation, no hyperparameter search, calibration with 299 patients is noisy.
