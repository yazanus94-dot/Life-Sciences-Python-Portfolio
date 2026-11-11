# Demo video script (4-5 min, screen recording with your voice)

Speak in your own words. Do not read this. Every section lists what is on screen.

**0:00 The question (20 s).** Screen: README title. "How does EEG frequency content change across sleep stages, and can I show it with statistics that hold up?"

**0:20 Data (30 s).** Screen: README data section. PhysioNet Sleep-EDF Expanded, Fpz-Cz channel, 100 Hz, hypnogram scored in 30-s epochs by human experts. Say why you picked one channel and the first night only.

**0:50 Method (90 s).** Screen: src/eeg_sleep. Walk through four files: preprocess.py (zero-phase 0.3-35 Hz filter, why zero phase; epoching; 500 uV peak-to-peak rejection and why not lower), spectral.py (Welch, 4-s Hann, relative power and why), analysis.py (subject-level medians, Friedman, paired Wilcoxon, Benjamini-Hochberg).

**2:20 Run it live (40 s).** Terminal: pytest, then python scripts/run_pipeline.py --subjects 0-9. Mention tests use a synthetic recording with a known injected frequency per stage, so correctness is checked without downloading data.

**3:00 Results (60 s).** Screen: hypnogram_spectrogram.png, relative_power_heatmap.png, stage_vs_wake_stats.csv. Describe ONLY what your run produced. Compare with textbook expectations (delta up in N3, spindle-range sigma up in N2, alpha highest in wake) and say plainly if anything differs.

**4:00 Evaluator problem (60 s).** Screen: problems/problem_01_pseudoreplication.md, then run problem_01_solution.py. "A model asked to compare N2 vs N3 will often pool epochs. In my simulation that gives a false-positive rate of 0.83 against a nominal 0.05; the subject-level test gives 0.04."

**5:00 Limitations (30 s).** Single channel, healthy-age cohort subset, first night only, amplitude threshold chosen by hand, no sensitivity analysis yet.

Tips: record at 1080p, zoom the terminal font, upload to YouTube as unlisted, put a linked thumbnail at the top of the README.
