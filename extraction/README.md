# Extraction scripts

Each script exports one source to the common specification in `../docs/data_specification.md`
and writes `<SOURCE>_creatinine_series.csv` and `<SOURCE>_episode_covariates.csv`, with
`<SOURCE>` exactly as in `analysis/config.R` (e.g. `MIMIC-IV`, `AmsterdamUMCdb`).
Edit the `path/to/...` placeholders at the top before running. Run each script whole
(`source()` or `Rscript`), not line by line: stepping through `data.table` code in RStudio
can trigger spurious errors.

| Source | Script | Written by | Status |
|---|---|---|---|
| MIMIC-IV | `MIMIC-IV_export.R` | core team | as run |
| INSPIRE | `INSPIRE_export.R` | core team | **changed**: t0 = operation start, not midnight of the operation date |
| NWICU | `NWICU_export.R` | core team, edited in use | as run |
| HiRID | `HiRID_export.R` | core team | as run, plus removal of repeated rows |
| SICdb | `SICdb_export.R` | research team | as run |
| Zigong | `Zigong_export.R` | core team | **corrected**: times are hours, not days (see below) |
| PIC | `PIC_export.R` | core team | as run (reads the V1.1.0 folder) |
| AmsterdamUMCdb | `AmsterdamUMCdb_export.R` | research team | as run |
| VitalDB | `VitalDB_export.R` | research team | rewritten to read the PhysioNet files; reproduces the published VitalDB results exactly |
| eICU-CRD | `eICU-CRD_export.R` | research team | as run (episode and patient counts match the paper); output names standardised, `unitvisitnumber` / `hospitaldischargeyear` added |
| MOVER (EPIC) | `MOVER_export.R` | research team | as run (episode and patient counts match the paper); comments translated, output names standardised, `t0_source` added |
| ZD, NJF, AHMU | `institutional/*_export.R` | research team | wrappers; **helper `export_creatinine.py` to be added** |

`institutional/template_institutional_export.R` is the specification template given to
the institutional data teams.

## Source-specific notes

- **HiRID** ships measurements as tar.gz archives; the script unpacks them once (about 10 GB).
  HiRID holds no pre-admission history and no discharge or death timestamp.
- **SICdb** records all times as offsets in seconds from ICU admission.
- **Zigong** records every time (LabTime, DISCHARGE_DATE_TIME, Death_Date) as an offset in
  **hours** from hospital admission (`datDictionary.csv`). The first version of the script
  read LabTime as days, which kept only the first 90 hours of results (7,613 of 16,031 in
  window) and made the "7-day" AKI window 7 hours. It also never linked `dtOutCome` (it
  joined on INP_NO, which dtOutCome does not carry), so no deaths were recorded.
- **INSPIRE** `opdate` is midnight of the operation day. Using it as t0 placed 9,915
  same-day pre-operative results after t0; t0 is now `opstart_time`.
- **HiRID** the published export repeated each draw 2-5 times (same second and value);
  the script and the pipeline now keep one row per draw.
- **VitalDB** creatinine is `name == "cr"` in `lab_data.csv`, already in mg/dL.
- **eICU-CRD** times are minutes from ICU admission (t0 = 0 for every stay). Results can be
  aligned only within one hospitalisation (`patienthealthsystemstayid`), so there is no
  pre-hospital history. Death = hospital discharge offset when `hospitaldischargestatus` is
  Expired.
- **MOVER** uses the EPIC subset only (64,357 operations, 39,685 patients). Datetimes are
  m/d/yy H:MM. t0 = entry to the operating room, else anaesthesia start, else SURGERY_DATE
  (midnight of the operation day) for 5,174 operations (8%); `t0_source` records which.
  In-hospital death is timed at discharge. `BIRTH_DATE` holds age in years.
- **PIC** ships two versions in one download. Read `V1.1.0/`; the root-level files are v1.0.0.
- **AmsterdamUMCdb** reads connection settings from environment variables (`AMDS_HOST`,
  `AMDS_PORT`, `AMDS_DB`, `AMDS_USER`). Set `PGPASSWORD` outside the script. Never write a
  password into this file.
- **NWICU** — derive in-hospital death from timestamps (death on or before discharge), not from
  `hospital_expire_flag`, which is coded differently from MIMIC-IV.

## Rules every export must follow

1. Export 365 days before to 90 days after t0, including outpatient and readmission results.
   Check the unit of every time offset against the source's data dictionary.
2. Declare the unit; convert umol/L to mg/dL (divide by 88.42) and say you did.
3. Export identifiers as text, never numeric.
4. Keep series and covariates as two files. Merging drops patients with no creatinine.
5. Never open the output in Excel.
