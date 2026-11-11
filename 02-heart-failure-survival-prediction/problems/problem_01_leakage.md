# Evaluator problem 1: target leakage through follow-up time

## Task (as given to a model)

Using the heart-failure clinical records (age, ejection fraction, serum creatinine,
..., follow-up "time", and DEATH_EVENT), build a classifier for DEATH_EVENT and
report its cross-validated AUC.

## Correct approach

1. Exclude "time" from the features. Follow-up ends early when a patient dies and
   runs to the end of the study otherwise, so time encodes the outcome. It is not
   available when the prediction would actually be made.
2. Better: model time-to-event directly with a Cox model and report a
   cross-validated concordance index, since follow-up lengths differ between patients.
3. Fit scalers and transforms inside each training fold (a Pipeline), not on the
   full dataset before splitting.
4. Report calibration and a baseline-referenced score (Brier vs prevalence), not AUC alone.

## Common model mistakes

| Mistake | Effect |
|---|---|
| Include time as a feature | AUC inflated by about 0.1 to 0.18 in the simulation |
| StandardScaler fit on all rows before CV | Mild leakage of test-fold statistics |
| Report AUC only | Hides miscalibration |
| Treat censored survivors as true negatives with no caveat | Conflates short follow-up with low risk |

## Verification

python problems/problem_01_solution.py (synthetic cohorts, same schema, event rate
about 38 percent, 5 seeds, repeated 5-fold CV). Output from my run:

- seed 0: 0.708 without time, 0.809 with (+0.101)
- seed 1: 0.719 / 0.863 (+0.143)
- seed 2: 0.735 / 0.830 (+0.095)
- seed 3: 0.674 / 0.798 (+0.124)
- seed 4: 0.665 / 0.846 (+0.181)

Run it on the real file with: python problems/problem_01_solution.py data/heart_failure_clinical_records_dataset.csv
and report what you get. The synthetic numbers are not findings about real patients.
