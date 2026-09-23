# Extraction scripts

Each script exports one source to the common specification in `../docs/data_specification.md`.
Edit the `path/to/...` placeholders at the top before running. Run each script whole
(`source()` or `Rscript`), not line by line: stepping through `data.table` code in RStudio
can trigger spurious errors.

| Source | Script | Written by | Status |
|---|---|---|---|
| MIMIC-IV | `MIMIC-IV_export.R` | core team | as run |
| INSPIRE | `INSPIRE_export.R` | core team | as run |
| NWICU | `NWICU_export.R` | core team, edited in use | as run |
| HiRID | `HiRID_export.R` | core team | as run |
| SICdb | `SICdb_export.R` | research team | as run |
| Zigong | `Zigong_export.R` | core team | as run |
| PIC | `PIC_export.R` | core team | as run (reads the V1.1.0 folder) |
| AmsterdamUMCdb | `AmsterdamUMCdb_export.R` | research team | as run |
| VitalDB | `VitalDB_export.R` | research team | as run |
| eICU-CRD | — | research team | **to be added** |
| MOVER | — | research team | **to be added** |
| ZD, NJF, AHMU | `institutional/*_export.R` | research team | wrappers; **helper `export_creatinine.py` to be added** |

`institutional/template_institutional_export.R` is the specification template given to
the institutional data teams.

## Source-specific notes

- **HiRID** ships measurements as tar.gz archives; the script unpacks them once (about 10 GB).
  HiRID holds no pre-admission history and no discharge or death timestamp.
- **SICdb** records all times as offsets in seconds from ICU admission.
- **Zigong** records lab times as offsets in days from hospital admission.
- **PIC** ships two versions in one download. Read `V1.1.0/`; the root-level files are v1.0.0.
- **AmsterdamUMCdb** reads connection settings from environment variables (`AMDS_HOST`,
  `AMDS_PORT`, `AMDS_DB`, `AMDS_USER`). Set `PGPASSWORD` outside the script. Never write a
  password into this file.
- **NWICU** — derive in-hospital death from timestamps (death on or before discharge), not from
  `hospital_expire_flag`, which is coded differently from MIMIC-IV.

## Rules every export must follow

1. Export 365 days before to 90 days after t0, including outpatient and readmission results.
2. Declare the unit; convert umol/L to mg/dL (divide by 88.42) and say you did.
3. Export identifiers as text, never numeric.
4. Keep series and covariates as two files. Merging drops patients with no creatinine.
5. Never open the output in Excel.
