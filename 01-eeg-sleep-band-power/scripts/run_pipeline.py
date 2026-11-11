"""Run the full analysis on Sleep-EDF subjects.

Usage: python scripts/run_pipeline.py --subjects 0-9 --out results
First run downloads about 50 MB per subject from PhysioNet (cached by MNE).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from eeg_sleep import analysis, plots  # noqa: E402
from eeg_sleep.config import STAGE_NAMES  # noqa: E402
from eeg_sleep.data import load_subject  # noqa: E402
from eeg_sleep.preprocess import bandpass, epoch_by_stage, reject_artifacts  # noqa: E402

SEED = 42


def parse_range(s: str) -> list[int]:
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default="0-9")
    ap.add_argument("--out", default="results")
    ap.add_argument("--channel", default="EEG Fpz-Cz")
    ap.add_argument("--ptp-uv", type=float, default=500.0)
    ap.add_argument("--synthetic", action="store_true", help="demo on synthetic recordings (no download)")
    a = ap.parse_args()
    np.random.seed(SEED)

    out = Path(a.out); (out / "figures").mkdir(parents=True, exist_ok=True)
    tables, qc, example = [], [], None
    for sid in parse_range(a.subjects):
        if a.synthetic:
            from eeg_sleep.synthetic import synthetic_raw
            raw = bandpass(synthetic_raw(sid))
        else:
            raw = bandpass(load_subject(sid, channel=a.channel))
        X, y = epoch_by_stage(raw)
        X, y, nrej = reject_artifacts(X, y, a.ptp_uv)
        sf = raw.info["sfreq"]
        tables.append(analysis.epoch_table(X, y, sf, sid))
        qc.append(dict(subject=sid, epochs_kept=len(y), epochs_rejected=nrej))
        if example is None:
            example = (X, y, sf)
        print(f"subject {sid}: kept {len(y)}, rejected {nrej}")

    epochs = pd.concat(tables, ignore_index=True)
    subj = analysis.subject_level(epochs)
    stats = analysis.stage_tests(subj)
    epochs.to_csv(out / "epoch_band_power.csv", index=False)
    subj.to_csv(out / "subject_stage_band_power.csv", index=False)
    stats.to_csv(out / "stage_vs_wake_stats.csv", index=False)
    pd.DataFrame(qc).to_csv(out / "qc_epoch_counts.csv", index=False)

    X, y, sf = example
    plots.hypnogram_and_spectrogram(X, y, sf, out / "figures" / "hypnogram_spectrogram.png")
    plots.mean_psd_by_stage(X, y, sf, out / "figures" / "psd_by_stage.png", STAGE_NAMES)
    plots.relative_power_heatmap(subj, out / "figures" / "relative_power_heatmap.png")
    plots.subject_boxplots(subj, out / "figures" / "subject_boxplots.png")
    print(stats.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
