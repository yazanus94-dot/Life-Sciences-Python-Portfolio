"""Shared constants. Bands are contiguous and half-open [lo, hi) so relative
power sums to exactly 1 across the analysed range."""

BANDS = {
    "delta": (0.5, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 12.0),
    "sigma": (12.0, 16.0),  # sleep-spindle range
    "beta": (16.0, 30.0),
}
TOTAL_RANGE = (0.5, 30.0)

# Sleep-EDF annotations -> integer stage codes. R&K stages 3 and 4 are merged
# into N3 following the AASM convention.
ANNOTATION_MAP = {
    "Sleep stage W": 0,
    "Sleep stage 1": 1,
    "Sleep stage 2": 2,
    "Sleep stage 3": 3,
    "Sleep stage 4": 3,
    "Sleep stage R": 4,
}
STAGE_NAMES = {0: "W", 1: "N1", 2: "N2", 3: "N3", 4: "REM"}
STAGE_ORDER = ["W", "N1", "N2", "N3", "REM"]

EPOCH_SEC = 30.0
