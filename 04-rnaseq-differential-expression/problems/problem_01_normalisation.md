# Evaluator problem 1: composition bias and total-count normalisation

## Task (as given to a model)

Two groups of RNA-seq samples. A handful of very highly expressed genes are
strongly up-regulated in the treated group. Write Python to normalise the counts
and test for differential expression. How many genes should be called?

## Correct approach

Use a composition-robust size factor (median-of-ratios, or TMM). Scaling by total
counts assumes the total RNA output is equal across samples. When a few dominant
genes absorb a larger share of reads in one group, every other gene looks
down-regulated in that group.

## Common model mistakes

| Mistake | Effect |
|---|---|
| Scale by library size / CPM, then test | Hundreds of spurious down-regulated genes |
| Run a t-test on raw counts | Ignores the mean-variance relationship of counts |
| Per-gene dispersion with 4 vs 4 samples, no shrinkage | Anti-conservative p-values (see model.py docstring) |
| Ignore the paired (subject) structure | Loses power from between-subject variance |

## Verification

python problems/problem_01_solution.py simulates 5 dominant genes up 8-fold in
treated samples and no other change (seed 7). Output from my run:

- total-count normalisation: 629 of 784 other genes falsely called (80.2 percent)
- median-of-ratios: 1 of 784 other genes falsely called (0.1 percent)
