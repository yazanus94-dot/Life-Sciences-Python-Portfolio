FEATURES = [
    "age", "anaemia", "creatinine_phosphokinase", "diabetes", "ejection_fraction",
    "high_blood_pressure", "platelets", "serum_creatinine", "serum_sodium", "sex", "smoking",
]
TIME, EVENT = "time", "DEATH_EVENT"
REQUIRED = FEATURES + [TIME, EVENT]
EF_CUTOFF = 40  # percent; clinical threshold for reduced ejection fraction (pre-specified)
SEED = 42
