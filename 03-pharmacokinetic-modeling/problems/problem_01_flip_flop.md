# Evaluator problem 1: flip-flop kinetics in oral PK

## Task (as given to a model)

Oral dose of 4 mg/kg; 11 concentration measurements over 24 h. Fit a
one-compartment model with first-order absorption (ka) and elimination (ke) and
report the elimination half-life ln(2)/ke.

## Correct approach

For oral data, the curve depends on ka and ke symmetrically: swapping ka and ke
and rescaling V/F to V/F * ke/ka gives exactly the same concentration-time curve.
The data alone cannot tell which rate is absorption and which is elimination.
Resolve it with outside information:

1. IV data (gives ke directly), or a known elimination half-life from the literature.
2. Pharmacological prior: for most drugs, absorption is faster than elimination.
3. The terminal log-linear slope equals the SMALLER of the two rates, so a fit with
   ke > ka means the apparent terminal phase is absorption-limited (flip-flop).
4. Always report both rates and flag the assumption.

## Common model mistakes

| Mistake | Effect |
|---|---|
| Report ln(2)/ke from whichever fit converged | Half-life wrong by more than an order of magnitude when the fit flips |
| Run a single start and trust it | Convergence depends on the start; both solutions have identical SSE |
| Call a good fit proof the model is right | Goodness of fit cannot distinguish the two solutions |
| Compare V/F across subjects without checking for flips | Flipped subjects have V/F off by ke/ka |

## Verification

python problems/problem_01_solution.py (seed 42, 5 percent multiplicative noise).
Output from my run:

- truth: ka 1.5, ke 0.1, V/F 0.45, half-life 6.93 h
- start with ka > ke: ka 1.583, ke 0.097, V/F 0.467, SSE 0.01990, half-life 7.12 h
- start with ka < ke: ka 0.097, ke 1.583, V/F 0.029, SSE 0.01990, "half-life" 0.44 h

Identical SSE, a 16-fold difference in half-life. tests/test_pk.py also proves the
swap identity analytically.
