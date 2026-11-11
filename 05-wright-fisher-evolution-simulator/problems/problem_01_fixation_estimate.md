# Evaluator problem 1: estimating a fixation probability by simulation

## Task (as given to a model)

A single new beneficial mutation (s = 0.01) arises in a diploid population of
N = 1000. Estimate by simulation the probability that it fixes, with a 95% confidence
interval, and compare with theory.

## Correct approach

1. Start from ONE mutant copy, p0 = 1/(2N), not from an intermediate frequency.
2. Use theory to size the experiment. Kimura gives u = 0.0198, so 100 replicates
   yield about 2 fixations. For a relative standard error near 2 percent you need
   thousands of fixations, i.e. about 10^5 replicates.
3. For a rare-event proportion use Wilson or Clopper-Pearson, not Wald.
4. Compare with Kimura (1 - exp(-4Nsp0)) / (1 - exp(-4Ns)), and note the 2s
   approximation is only accurate for small s and large N.

## Common model mistakes

| Mistake | Effect |
|---|---|
| 100 replicates and a Wald interval | 13.5 percent chance of zero fixations; Wald then returns [0, 0]; coverage is below nominal |
| Use 1/(2N) as the answer for a selected allele | 40-fold too low (0.0005 vs 0.0198) |
| Start at p = 0.5 | Measures something else; no longer a new mutation |
| Treat "inside the CI" as proof when the CI is huge | Uninformative; a wide interval also contains wrong theories |

## Verification

python problems/problem_01_solution.py (seed 42). Output from my run:

- Kimura u = 0.01980 (2s = 0.020; neutral 0.0005)
- coverage of nominal 95% intervals with n = 100: Wald 0.862, Wilson 0.951, Clopper-Pearson 0.984
- P(zero fixations in 100 replicates) = 0.135
- 200000 replicates: 3899 fixations, estimate 0.01949, 95% CI [0.01889, 0.02011]; Kimura inside the CI
