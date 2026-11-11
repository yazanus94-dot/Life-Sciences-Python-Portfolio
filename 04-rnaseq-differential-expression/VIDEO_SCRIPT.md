# Demo video script (4-5 min, screen recording with your voice)

Speak in your own words. Do not read this.

**0:00 Question (20 s).** Which genes respond to dexamethasone in human airway smooth muscle cells, and does my pipeline recover the known biology?

**0:20 Data (30 s).** Himes et al. 2014, GEO GSE52778, Bioconductor airway package: 8 samples, 4 donors, each untreated and dexamethasone-treated. Explain why this is a paired design and why that matters.

**0:50 Method (100 s).** Walk through src/rnaseq_de: filter_low_counts; median_of_ratios (why not total counts); the negative-binomial GLM with donor + treatment; dispersion shrinkage and why per-gene estimates are unreliable with 8 samples; Benjamini-Hochberg. Be honest that this is a simplified DESeq2.

**2:30 Run it live (30 s).** pytest (simulated truth: recall, FDR, null behaviour), then run_pipeline.py on the airway files.

**3:00 Results (60 s).** Show volcano, PCA (does donor or treatment separate?), known_gene_validation.csv. Report ONLY what your run produced. If a known gene is not recovered, say so and offer a hypothesis.

**4:00 Evaluator problem (45 s).** Run problem_01_solution.py: total-count normalisation falsely calls 80 percent of genes in that simulation; median-of-ratios calls 0.1 percent.

**4:45 Limitations (30 s).** Simplified dispersion model, not cross-checked vs PyDESeq2 until you do it, one dataset, no pathway analysis, gene-symbol collapsing.
