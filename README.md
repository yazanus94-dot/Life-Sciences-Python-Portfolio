# Life Sciences Python Portfolio

A curated portfolio of five end-to-end Python projects spanning biomedical signal analysis, clinical survival/risk modeling, pharmacokinetics, RNA-seq differential expression, and population-genetics simulation. Each project includes source code, tests, reproducible scripts, documentation, and evaluation/interpretation artifacts designed to demonstrate rigorous, evidence-oriented scientific computing.

**Author:** [Yazan Aljabari](https://www.linkedin.com/in/yazan-aljabari-878202139/)  
**License:** [MIT](LICENSE)  
**Branch:** `main`

---

## Overview

This repository brings together applied life-sciences workflows that emphasize correct statistical units of inference, leakage-safe modeling, model validation against analytic/theoretical expectations, and transparent reporting of limitations. The projects are structured as independent, reproducible modules that can be explored, tested, and extended individually.

### Included projects

| # | Project | Domain | Core question / goal |
|---|---|---|---|
| 1 | [EEG Sleep Band Power](01-eeg-sleep-band-power/) | Neurophysiology / signal processing | How does relative EEG band power differ across sleep stages when the subject (not epoch) is the unit of inference? |
| 2 | [Heart Failure Survival Prediction](02-heart-failure-survival-prediction/) | Clinical biostatistics / ML | Which clinical factors are associated with death after heart failure, and how well can risk be predicted without target leakage? |
| 3 | [Pharmacokinetic Modeling](03-pharmacokinetic-modeling/) | Pharmacokinetics | What are oral theophylline absorption/elimination kinetics and inter-subject variability, and which compartment model is supported? |
| 4 | [RNA-seq Differential Expression](04-rnaseq-differential-expression/) | Transcriptomics | Which genes respond to dexamethasone in a paired airway design, and does the pipeline recover known biology? |
| 5 | [Wright-Fisher Evolution Simulator](05-wright-fisher-evolution-simulator/) | Population genetics | Does simulation reproduce analytic fixation, fixation-time, and heterozygosity-decay predictions under drift and selection? |

---

## Problem statement

Life-sciences data analysis is frequently undermined by subtle but critical errors: pseudoreplication, data leakage, inappropriate normalization, weak model identification, and overconfident interpretation of small or simulated samples. This portfolio addresses that gap by implementing compact, inspectable pipelines that pair domain methods with explicit correctness checks (evaluator problems, statistical validation, and reproducibility controls).

## Objectives

- Demonstrate practical, end-to-end Python workflows across multiple life-sciences domains.
- Apply statistically defensible methodology (correct inference units, FDR control, leakage-safe CV, calibration/validation, theory-consistent simulation).
- Make results interpretable and auditable through documented assumptions, included artifacts, and tests.
- Provide a professional, recruiter- and reviewer-friendly repository structure with clear setup, usage, and limitations.

## Features

- **Five self-contained projects** with `src/`, `scripts/`, `tests/`, `docs/`, and project-specific documentation.
- **Methodological rigor:** subject-level inference (EEG), pipeline-contained preprocessing/CV (heart failure), AICc/NCA cross-checks (PK), median-of-ratios normalization & FDR (RNA-seq), Clopper–Pearson/Wilson intervals and analytic validation (Wright–Fisher).
- **Evaluator/problem modules** in each project illustrating common failure modes (pseudoreplication, leakage, flip-flop kinetics, composition bias, fixation-estimate uncertainty).
- **Reproducibility support:** pinned/declared dependencies, scripted pipelines, synthetic/no-download demo modes where applicable, and automated tests.
- **Included artifacts:** real-data CSV results and figures for pharmacokinetic and Wright–Fisher runs included in-repo; synthetic/demo media clearly labeled where used.
- **CI workflow** for repository-level validation.

## Technologies

- **Language:** Python (3.10+ recommended)
- **Scientific stack:** NumPy, SciPy, pandas, matplotlib, seaborn, statsmodels
- **Domain/ML tools:** scikit-learn, lifelines, MNE-Python
- **Testing/quality:** pytest, GitHub Actions
- **Data formats:** CSV, EDF/annotations (via MNE when downloaded), figures/media (PNG/GIF/MP4)

## Project structure

```text
.
├── 01-eeg-sleep-band-power/                # Sleep-EDF spectral analysis (MNE), subject-level stats/FDR
├── 02-heart-failure-survival-prediction/   # KM/Cox + leakage-safe classification/calibration
├── 03-pharmacokinetic-modeling/            # 1/2-cmt oral PK fitting, AICc, NCA, bootstrap, superposition
├── 04-rnaseq-differential-expression/      # Median-of-ratios, NB-GLM, paired design, FDR, validation
├── 05-wright-fisher-evolution-simulator/   # Drift/selection simulation validated against Kimura theory
├── docs/                                   # Portfolio-level results/dashboard assets
├── .github/workflows/ci.yml                # Continuous integration
├── index.html                              # Portfolio landing page
├── requirements.txt                        # Consolidated Python dependencies
├── pyproject.toml                          # Project/test configuration
├── LICENSE                                 # MIT license
└── README.md                               # This file
```

Each project directory typically contains: `README.md`, `src/`, `scripts/`, `tests/`, `docs/` (media), `problems/` (evaluator scenarios), and project-specific config/data/results as applicable.

## Installation

### Prerequisites
- Python **3.10+**
- `pip` and virtual environment tooling
- (Optional, project-specific) R + Bioconductor for `scripts/export_airway.R` in the RNA-seq project; `ffmpeg`/`pillow` to re-render demo videos; internet access for first-run dataset downloads (Sleep-EDF via MNE, UCI heart-failure CSV if not provided locally)

### Setup

```bash
# Clone
git clone https://github.com/yazanus94-dot/Life-Sciences-Python-Portfolio.git
cd Life-Sciences-Python-Portfolio

# Create and activate a virtual environment
python -m venv .venv
# Windows (PowerShell): .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

> **Note:** Some projects ship their own `requirements.txt`/`pyproject.toml` for isolated use; the root `requirements.txt` provides a consolidated environment covering the portfolio.

## Usage

Run commands from the repository root (or `cd` into a project) with `PYTHONPATH`/editable context as each project expects (most scripts assume execution from the project directory or include `src` on path). Examples:

```bash
# 1) EEG sleep band power (synthetic/no-download demo)
cd 01-eeg-sleep-band-power
python scripts/run_pipeline.py --synthetic --subjects 0-7

# 2) Heart failure (synthetic demo cohort)
cd ../02-heart-failure-survival-prediction
python scripts/run_pipeline.py --synthetic

# 3) Pharmacokinetic modeling (uses included data/Theoph.csv)
cd ../03-pharmacokinetic-modeling
python scripts/run_pipeline.py

# 4) RNA-seq differential expression (synthetic demo; no R required)
cd ../04-rnaseq-differential-expression
python scripts/make_synthetic_data.py
python scripts/run_pipeline.py --counts data/synthetic_counts.csv --coldata data/synthetic_coldata.csv --factors subject treat --contrast treat --ref ctl --skip-validation

# 5) Wright–Fisher evolution simulator (pure simulation)
cd ../05-wright-fisher-evolution-simulator
python scripts/run_pipeline.py
```

For real-data runs, follow each project README (e.g., Sleep-EDF download via MNE, UCI heart-failure CSV placement, airway count export via R/Bioconductor). Synthetic/demo outputs are **not** real biological/clinical findings and are labeled as such in project docs.

## Methodology

| Project | Data / inputs | Core methods | Validation & safeguards |
|---|---|---|---|
| EEG Sleep Band Power | PhysioNet Sleep-EDF Expanded (downloaded via MNE) or synthetic recordings | FIR band-pass, 30-s epoching by hypnogram, Welch PSD, relative band power, subject-level medians | Friedman/Wilcoxon with BH-FDR; pseudoreplication evaluator; artifact rejection; offline unit tests (sinusoids, Parseval, known injections) |
| Heart Failure Survival Prediction | UCI Heart Failure Clinical Records (299 pts) or synthetic cohort | Schema validation, KM/log-rank, penalized Cox (PH checks), repeated stratified CV (Cox C-index, logreg/RF AUC), calibration/Brier | Transforms/scaling inside CV pipelines; explicit leakage demo (follow-up time); schema/sign/effect tests |
| Pharmacokinetic Modeling | Theophylline dataset (`data/Theoph.csv`, included) | 1-/2-compartment oral models (closed-form & ODE), multi-start log-scale NLS, AICc, CL/F/AUC/t½, NCA cross-check, subject bootstrap, superposition vs steady-state | ODE vs analytic agreement tests; identifiability/flip-flop evaluator; included per-subject and population summary artifacts |
| RNA-seq Differential Expression | Airway (GSE52778) via R export or synthetic paired counts | Low-count filtering, median-of-ratios size factors, NB-GLM with paired design, dispersion shrinkage, Wald tests + BH-FDR, PCA/plots | Normalization-bias evaluator; recall/FDR/null calibration tests on synthetic truth; known-gene validation when real data available |
| Wright–Fisher Evolution Simulator | In-silico (no external data) | Vectorized WF with drift/selection, variable N (bottlenecks), absorbing states | Statistical agreement with Kimura fixation/time and heterozygosity theory; Wald/Wilson/Clopper–Pearson interval comparisons; reproducibility tests |

## Results

> **Transparency:** This section reports only artifacts included in this repository or explicitly labeled demo/simulation outputs. It does **not** fabricate real experimental results. Projects requiring external downloads (EEG real Sleep-EDF, UCI real CSV, airway real counts) produce user-run results under `results/` as documented in each project README.

### Included run artifacts (from repository files)

#### Pharmacokinetic Modeling — population summary (`03-pharmacokinetic-modeling/results/population_summary.csv`)
Geometric means with subject-level bootstrap 95% CI (included run on `data/Theoph.csv`):

| Parameter | Geometric mean | 95% CI (low) | 95% CI (high) |
|---|---:|---:|---:|
| ka (1/h) | 1.359 | 0.898 | 2.230 |
| ke (1/h) | 0.088 | 0.078 | 0.097 |
| V/F (L/kg) | 0.454 | 0.419 | 0.491 |
| CL/F (L/h/kg) | 0.040 | 0.034 | 0.046 |
| t½ (h) | 7.898 | 7.165 | 8.904 |

Additional included artifacts: `results/subject_fits.csv`, `results/figures/subject_fits.png`, `results/figures/accumulation.png`.

#### Wright–Fisher Evolution Simulator — fixation probabilities (`05-wright-fisher-evolution-simulator/results/fixation_probability.csv`)
Included simulation (N = 50, 20,000 replicates per *s*, seed 42):

| s | Fixations / reps | Estimated P(fix) | 95% CI (Clopper–Pearson) | Kimura theory | Theory in CI | Mean t(fix) |
|---:|---:|---:|---:|---:|:---:|---:|
| 0.00 | 222 / 20000 | 0.0111 | [0.0097, 0.0127] | 0.0100 | Yes | 206.30 |
| 0.01 | 462 / 20000 | 0.0231 | [0.0211, 0.0253] | 0.0229 | Yes | 195.65 |
| 0.02 | 791 / 20000 | 0.0396 | [0.0369, 0.0423] | 0.0399 | Yes | 162.48 |
| 0.05 | 1908 / 20000 | 0.0954 | [0.0914, 0.0996] | 0.0952 | Yes | 108.75 |
| 0.10 | 3532 / 20000 | 0.1766 | [0.1713, 0.1820] | 0.1813 | Yes | 71.30 |

Additional included artifacts: `results/figures/fixation_vs_s.png`, `results/figures/heterozygosity.png`, `results/figures/trajectories.png`.

#### Other projects
- **EEG Sleep Band Power / Heart Failure / RNA-seq:** repositories include pipelines, tests, evaluator problems, and demo media; quantitative findings for real datasets are produced by running the pipelines on the respective real inputs (see each project README). Any `--synthetic` runs or demo videos generated from synthetic data are explicitly **non-clinical/non-biological demos**.
- See [docs/results.md](docs/results.md) for a consolidated dashboard view of included artifacts.

## Testing

Automated tests are included per project (offline/synthetic where possible).

```bash
# Run all tests from repository root (uses pyproject.toml testpaths)
pytest

# Or run per project
cd 01-eeg-sleep-band-power && pytest
cd 02-heart-failure-survival-prediction && pytest
cd 03-pharmacokinetic-modeling && pytest
cd 04-rnaseq-differential-expression && pytest
cd 05-wright-fisher-evolution-simulator && pytest
```

CI executes the test suite via GitHub Actions (`.github/workflows/ci.yml`).

## Limitations

- **Data access/licensing:** Some real datasets require download/acceptance of external terms (PhysioNet Sleep-EDF, UCI heart-failure CSV, GEO/Bioconductor airway). The repository does not redistribute restricted raw data.
- **Generalizability:** Included cohorts/datasets are small or domain-specific (e.g., 12 PK subjects, 299 HF patients, 8 RNA-seq samples), limiting statistical power and external validity.
- **Modeling simplifications:** Per-subject (non-NLME) PK fits; binary HF classification ignores censoring (Cox is primary); simplified NB dispersion/shrinkage vs DESeq2; single-locus WF model without mutation/recombination/linkage/structure.
- **Reproducibility environment:** Numeric results can vary slightly across platforms/library versions; demo/synthetic outputs are illustrative only.
- **No clinical/diagnostic use:** This portfolio is for educational and methodological demonstration, not medical decision-making.

## Future improvements

- Add containerized environments (Docker) and pinned lockfiles for fully deterministic reruns.
- Integrate richer reporting notebooks and automated result manifests for user-executed real-data runs.
- Expand cross-validation/external validation where additional public datasets are available.
- Benchmark RNA-seq dispersion/Wald implementation against PyDESeq2/DESeq2 outputs on reference data.
- Extend PK analysis to population (NLME) modeling and covariate effects; add sensitivity analyses for EEG artifact thresholds.
- Enhance CI with linting/type checks and selective heavy-dependency job matrices.

## References

- Kemp B, et al. *Analysis of a sleep-dependent neuronal feedback loop: the slow-wave microcontinuity of the EEG.* IEEE Trans Biomed Eng, 2000. (Sleep-EDF)
- Goldberger AL, et al. *PhysioBank, PhysioToolkit, and PhysioNet.* Circulation, 2000.
- Ahmad T, et al. *Survival analysis of heart failure patients: A case study.* PLoS ONE, 2017. (UCI Heart Failure Clinical Records)
- Chicco D, Jurman G. *Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone.* BMC Med Inform Decis Mak, 2020.
- Boeckmann AJ, Sheiner LB, Beal SL. *NONMEM Users Guide*, 1994. (Theophylline dataset as distributed in R datasets/Rdatasets)
- Himes BE, et al. *RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene...* PLoS ONE, 2014. (GEO GSE52778; Bioconductor airway)
- Kimura M. *On the probability of fixation of mutant genes in a population.* Genetics, 1962.
- Crow JF, Kimura M. *An Introduction to Population Genetics Theory.* 1970.

## License

This repository is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

**Copyright (c) 2026 Yazan Aljabari** — [LinkedIn](https://www.linkedin.com/in/yazan-aljabari-878202139/)