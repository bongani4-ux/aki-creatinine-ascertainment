# Common export specification

## `<SOURCE>_creatinine_series.csv` — one row per creatinine result

| Column | Type | Description |
|---|---|---|
| `episode_id` | text | Joins to the covariate file |
| `patient_id` | text | Identifies repeat episodes |
| `t_days` | numeric | Days from t0; negative = before the index episode |
| `value` | numeric | Serum creatinine, mg/dL (the pipeline stops if `unit` is anything else) |
| `unit` | text | Always `mg/dL` after conversion |
| `charttime` / `timestamp` | datetime | Clock time of the draw, where the source records it |

Window: t_days from -365 to +90. No value is pre-selected as a baseline.
Exact duplicate rows (same episode, t_days and value) are removed by the pipeline
(`DEDUPLICATE` in `analysis/config.R`).

## `<SOURCE>_episode_covariates.csv` — one row per episode

| Column | Description |
|---|---|
| `episode_id`, `patient_id` | as above |
| `t0` (or `t0_min`) | index time: ICU admission (critical care; hospital admission for Zigong) or operation start (perioperative, surgical). A datetime, or a numeric offset on the source's own time axis |
| `age`, `sex` | demographics |
| `discharge_t_days` or `discharge_time` | discharge, as days from t0 or as a datetime |
| `death_t_days` or `death_time` | death, as days from t0 or as a datetime; empty if alive |

Absent fields are left present and empty, never dropped.
