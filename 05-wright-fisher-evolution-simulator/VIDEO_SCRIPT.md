# Demo video script (4-5 min, screen recording with your voice)

Speak in your own words. Do not read this. You can show docs/demo.gif on screen.

**0:00 Question (20 s).** How strongly do drift and selection shape the fate of a new mutation, and can a simulation reproduce population-genetics theory?

**0:20 Model (60 s).** src/wfsim/sim.py: 2N gene copies, selection step p' = p(1+s)/(1+sp), then Binomial(2N, p') drift; absorbing at 0 and 1. Explain why the active-set loop makes 10^5 replicates cheap.

**1:20 Theory (40 s).** analytic.py: neutral fixation probability 1/(2N); Kimura's formula and the 2s limit; conditional fixation time 4N; heterozygosity decay and the harmonic mean of N.

**2:00 Run it live (40 s).** pytest (statistical tests use a 4-standard-error tolerance, not arbitrary thresholds), then python scripts/run_pipeline.py.

**2:40 Results (70 s).** Show fixation_vs_s.png and the CSV. Report ONLY what your run produced: which s values have the Kimura prediction inside the CI; note the borderline case at large s and why (the formula is a diffusion approximation). Show the bottleneck heterozygosity figure.

**3:50 Evaluator problem (60 s).** Run problem_01_solution.py: 100 replicates and a Wald interval cover the truth 86 percent of the time and can return [0, 0]; Wilson 95 percent; 200000 replicates give a tight estimate.

**4:50 Limitations (20 s).** Single locus, no mutation, no linkage, no sex or age structure, haploid-style selection (genic).
