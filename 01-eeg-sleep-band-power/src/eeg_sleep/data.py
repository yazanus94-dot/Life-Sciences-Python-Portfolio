"""Loading Sleep-EDF (PhysioNet) via MNE downloader."""
from __future__ import annotations

import mne
from mne.datasets.sleep_physionet.age import fetch_data


def load_subject(
    subject: int,
    recording: int = 1,
    channel: str = "EEG Fpz-Cz",
    wake_buffer_min: float = 30.0,
) -> mne.io.BaseRaw:
    """Download (cached) and load one Sleep-EDF night as a single-channel Raw.

    Wake periods before sleep onset and after final awakening are cropped to
    wake_buffer_min minutes so wake does not dominate the epoch counts
    (same convention as the MNE sleep-staging tutorial).
    """
    psg_file, hyp_file = fetch_data(subjects=[subject], recording=[recording])[0]
    raw = mne.io.read_raw_edf(psg_file, preload=True, verbose="ERROR")
    ann = mne.read_annotations(hyp_file)
    buf = wake_buffer_min * 60.0
    if len(ann) > 2:
        ann.crop(ann[1]["onset"] - buf, ann[-2]["onset"] + buf)
    raw.set_annotations(ann, emit_warning=False)
    raw.pick([channel])
    return raw
