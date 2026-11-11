# Evaluator problem 1: pseudoreplication in sleep EEG statistics

## Task (as given to a model)

You have 30-s EEG epochs from 10 subjects (about 200 epochs per subject per
stage), with relative delta power computed per epoch. Test whether delta power
differs between N2 and N3 sleep. Write Python that returns a p-value.

## Correct approach

1. Aggregate to one value per subject per stage (median of that subject epochs).
2. Run a paired test across subjects (Wilcoxon signed-rank, or paired t-test if
   differences look normal). The subject is the unit of inference.
3. If testing several bands or stages, control FDR (Benjamini-Hochberg) across
   the whole family of tests.

## Common model mistakes

| Mistake | Why it is wrong |
|---|---|
| Pool all epochs and run an unpaired test (2000 vs 2000 "samples") | Epochs within a night are autocorrelated and share the subject baseline. Effective sample size is about 10, not 2000. |
| Test each band separately and report raw p-values | Inflates false discoveries across 5 bands x 4 stages. |
| Use absolute power for group comparison | Dominated by between-subject amplitude differences (skull thickness, electrode impedance). Relative power removes this. |

## Verification

Run python problems/problem_01_solution.py. It simulates a null experiment (no
true average stage effect, but each subject has their own random difference)
with seed 42 and 2000 simulations. Output from my run:

- epoch-level pooled test false-positive rate: 0.831
- subject-level paired test false-positive rate: 0.043
- nominal alpha: 0.05

The pooled test declares a "significant" difference in most null experiments.
