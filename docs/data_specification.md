# Common export specification

## `<SOURCE>_creatinine_series.csv` — one row per creatinine result

| Column | Type | Description |
|---|---|---|
| `episode_id` | text | Joins to the covariate file |
| `patient_id` | text | Identifies repeat episodes |
| `t_days` | numeric | Days from t0; negative = before the index episode |
| `value` | numeric | Serum creatinine, mg/dL |
| `unit` | text | Always `mg/dL` after conversion |
| `charttime` / `timestamp` | datetime | Clock time of the draw, where the source records it |

Window: t_days from -365 to +90. No value is pre-selected as a baseline.

## `<SOURCE>_episode_covariates.csv` — one row per episode

| Column | Description |
|---|---|
| `episode_id`, `patient_id` | as above |
| `t0` | index time: ICU admission (critical care) or operation (perioperative, surgical) |
| `age`, `sex` | demographics |
| `discharge_t_days` or `discharge_time` | discharge, as days from t0 or as a datetime |
| `death_t_days` or `death_time` | death, as days from t0 or as a datetime; empty if alive |

Absent fields are left present and empty, never dropped.
