"""Reference solution for problems/problem_01_pseudoreplication.md.

Simulates a null experiment (NO true stage effect on average) in which each
subject has their own random condition difference. Compares:
  A) epoch-level test: all epochs pooled as if independent  (wrong)
  B) subject-level paired test on per-subject medians       (correct)
and reports the false-positive rate at alpha = 0.05.
"""
import numpy as np
from scipy.stats import ttest_ind, ttest_rel

rng = np.random.default_rng(42)
N_SIM, N_SUBJ, N_EPOCH, SUBJ_SD, EPOCH_SD, ALPHA = 2000, 10, 200, 1.0, 1.0, 0.05

fp_epoch = fp_subj = 0
for _ in range(N_SIM):
    mu_a = rng.normal(0, SUBJ_SD, N_SUBJ)  # subject-specific level in stage A
    mu_b = rng.normal(0, SUBJ_SD, N_SUBJ)  # subject-specific level in stage B
    a = mu_a[:, None] + rng.normal(0, EPOCH_SD, (N_SUBJ, N_EPOCH))
    b = mu_b[:, None] + rng.normal(0, EPOCH_SD, (N_SUBJ, N_EPOCH))
    fp_epoch += ttest_ind(a.ravel(), b.ravel()).pvalue < ALPHA
    fp_subj += ttest_rel(np.median(a, 1), np.median(b, 1)).pvalue < ALPHA

print(f"False-positive rate, epoch-level pooled : {fp_epoch / N_SIM:.3f}")
print(f"False-positive rate, subject-level paired: {fp_subj / N_SIM:.3f}")
print(f"Nominal alpha                            : {ALPHA}")
