rm(list = ls())

library(data.table)

# -----------------------------------------------------------------------------
# SICdb creatinine-series export
#
# SICdb is pseudonymized and provides relative offsets rather than calendar
# datetimes. In both output files, timestamp-like fields are therefore numeric
# seconds from the patient's first primary MetaVision admission:
#   timestamp, t0, discharge_time, death_time
# t0 is the first transfer to an ICU/intermediate-care unit for that CaseID.
#
# Cross-readmission coverage is retained by aligning every CaseID on
# OffsetAfterFirstAdmission and then matching all of a patient's creatinine
# results to every eligible episode in the requested -365 to +90 day window.
# -----------------------------------------------------------------------------

sicdb_dir <- "path/to/SICdb"
# Set SICDB_OUT_DIR before running to override this destination. An ASCII-only
# default avoids path-encoding failures in older Windows R installations.
out_dir   <- Sys.getenv("SICDB_OUT_DIR", unset = "path/to/sicdb_output")

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90
SEC_PER_DAY   <- 86400

# Serum/plasma creatinine assays in the SICdb reference dictionary.
# 339 currently has no rows in laboratory.csv but is retained for compatibility.
# Excluded on purpose: urine (340, 370), HPM urine-panel (369), punctate fluid
# (380), clearance (441), and calcium/creatinine ratio (582).
SERUM_CREATININE_IDS <- c(339L, 367L, 368L)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

resolve <- function(dir, stem) {
  candidates <- c(
    file.path(dir, paste0(stem, ".csv")),
    file.path(dir, paste0(stem, ".csv.gz")),
    file.path(dir, stem, paste0(stem, ".csv"))
  )
  for (path in candidates) {
    if (file.exists(path) && !dir.exists(path)) {
      cat(sprintf("  %-14s -> %s\n", stem, path))
      return(path)
    }
  }
  stop(
    "Could not find a readable file for '", stem, "' under:\n  ", dir,
    "\nTried:\n  ", paste(candidates, collapse = "\n  ")
  )
}

require_columns <- function(dt, required, table_name) {
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(
      table_name, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

cat("=== Resolving SICdb input files ===\n")
f_cases <- resolve(sicdb_dir, "cases")
f_labs  <- resolve(sicdb_dir, "laboratory")
f_refs  <- resolve(sicdb_dir, "d_references")

# -----------------------------------------------------------------------------
# 1. Episode-level covariates
cat("\n=== Building episode covariates ===\n")

cases <- fread(f_cases)
refs  <- fread(f_refs)

require_columns(
  cases,
  c(
    "CaseID", "PatientID", "AdmissionYear", "TimeOfStay", "ICUOffset",
    "OffsetOfDeath", "Sex", "OffsetAfterFirstAdmission"
  ),
  "cases"
)
require_columns(
  refs,
  c("ReferenceGlobalID", "ReferenceValue", "ReferenceName", "ReferenceUnit"),
  "d_references"
)

if (anyDuplicated(cases$CaseID)) stop("CaseID is not unique in cases.csv.")

decode_reference <- function(x) {
  as.character(refs$ReferenceValue[match(as.integer(x), refs$ReferenceGlobalID)])
}

case_origin_sec <- suppressWarnings(as.numeric(cases$OffsetAfterFirstAdmission))
icu_offset_sec  <- suppressWarnings(as.numeric(cases$ICUOffset))
stay_sec        <- suppressWarnings(as.numeric(cases$TimeOfStay))
death_offset    <- suppressWarnings(as.numeric(cases$OffsetOfDeath))

invalid_icu_los <- is.na(stay_sec) | is.na(icu_offset_sec) | stay_sec < icu_offset_sec

episodes <- data.table(
  source          = "SICdb",
  episode_id      = as.character(cases$CaseID),
  patient_id      = as.character(cases$PatientID),
  t0              = case_origin_sec + icu_offset_sec,
  age             = suppressWarnings(as.numeric(cases$AgeOnAdmission)),
  sex             = decode_reference(cases$Sex),
  asa             = NA_real_,
  emergency       = NA_integer_,
  hosp_los_days   = stay_sec / SEC_PER_DAY,
  icu_los_days    = (stay_sec - icu_offset_sec) / SEC_PER_DAY,
  discharge_time  = case_origin_sec + stay_sec,
  death_time      = case_origin_sec + death_offset
)

episodes[invalid_icu_los, icu_los_days := NA_real_]
episodes[age < 0 | age > 110, age := NA_real_]

if (anyNA(episodes$t0)) stop("At least one episode has no usable t0.")
if (anyDuplicated(episodes$episode_id)) stop("episode_id is not unique.")

covariate_path <- file.path(out_dir, "SICdb_episode_covariates.csv")
fwrite(episodes, covariate_path, na = "")
cat(sprintf("  episodes written: %s\n", format(nrow(episodes), big.mark = ",")))
cat("  ASA: unavailable in the public SICdb tables; exported as missing.\n")
cat("  emergency status: no equivalent public field; exported as missing.\n")

# -----------------------------------------------------------------------------
# 2. Complete serum-creatinine series across all linked admissions
cat("\n=== Extracting serum creatinine ===\n")

labs <- fread(
  f_labs,
  select = c("id", "CaseID", "LaboratoryID", "Offset", "LaboratoryValue")
)
require_columns(
  labs,
  c("id", "CaseID", "LaboratoryID", "Offset", "LaboratoryValue"),
  "laboratory"
)

creat_ref <- refs[
  ReferenceGlobalID %in% SERUM_CREATININE_IDS,
  .(
    LaboratoryID = as.integer(ReferenceGlobalID),
    laboratory_name = as.character(ReferenceValue),
    raw_unit = as.character(ReferenceUnit)
  )
]
if (!nrow(creat_ref)) stop("No selected creatinine IDs exist in d_references.csv.")

cat("  selected reference entries:\n")
print(creat_ref)

cr <- labs[
  LaboratoryID %in% SERUM_CREATININE_IDS &
    is.finite(LaboratoryValue) & LaboratoryValue > 0,
  .(
    lab_row_id = as.character(id),
    lab_case_id = as.character(CaseID),
    LaboratoryID = as.integer(LaboratoryID),
    lab_offset_sec = suppressWarnings(as.numeric(Offset)),
    raw_value = suppressWarnings(as.numeric(LaboratoryValue))
  )
]

cr <- merge(cr, creat_ref, by = "LaboratoryID", all.x = TRUE, sort = FALSE)
if (anyNA(cr$raw_unit)) stop("A selected creatinine assay has no unit mapping.")

unit_key <- tolower(trimws(cr$raw_unit))
is_umol  <- grepl("umol|µmol|μmol", unit_key)
is_mgdl  <- grepl("mg\\s*/?\\s*dl", unit_key)
if (any(!is_umol & !is_mgdl)) {
  stop(
    "Unsupported creatinine unit(s): ",
    paste(unique(cr$raw_unit[!is_umol & !is_mgdl]), collapse = ", ")
  )
}

# Creatinine conversion: mg/dL = micromol/L / 88.4.
cr[, value := fifelse(is_umol, raw_value / 88.4, raw_value)]

case_axis <- data.table(
  lab_case_id = as.character(cases$CaseID),
  patient_id = as.character(cases$PatientID),
  case_origin_sec = case_origin_sec
)

unmatched_lab_cases <- setdiff(unique(cr$lab_case_id), case_axis$lab_case_id)
if (length(unmatched_lab_cases)) {
  stop(length(unmatched_lab_cases), " creatinine CaseID values do not match cases.csv.")
}

cr <- merge(cr, case_axis, by = "lab_case_id", all.x = TRUE, sort = FALSE)
cr[, timestamp := case_origin_sec + lab_offset_sec]

# Every result is matched to every ICU episode of the same patient. The absolute
# relative time axis permits inclusion of preceding/subsequent SICdb admissions.
ser <- merge(
  cr[, .(
    patient_id, lab_row_id, lab_case_id, timestamp,
    value, converted_from_umol = is_umol
  )],
  episodes[, .(episode_id, patient_id, t0)],
  by = "patient_id",
  allow.cartesian = TRUE,
  sort = FALSE
)
ser[, t_days_unrounded := (timestamp - t0) / SEC_PER_DAY]
ser <- ser[
  t_days_unrounded >= -LOOKBACK_DAYS &
    t_days_unrounded <= FOLLOWUP_DAYS
]

if (anyDuplicated(ser[, .(episode_id, lab_row_id)])) {
  stop("A laboratory row was duplicated within an episode after joining.")
}

series <- ser[, .(
  source = "SICdb",
  episode_id,
  patient_id,
  timestamp,
  t_days = round(t_days_unrounded, 5),
  value = round(value, 6),
  unit = "mg/dL"
)]
setorder(series, episode_id, t_days, timestamp)

series_path <- file.path(out_dir, "SICdb_creatinine_series.csv")
fwrite(series, series_path, na = "")
cat(sprintf("  creatinine rows written: %s\n", format(nrow(series), big.mark = ",")))

# -----------------------------------------------------------------------------
# 3. Validation
cat("\n=== Validation ===\n")

if (!identical(names(series), c(
  "source", "episode_id", "patient_id", "timestamp", "t_days", "value", "unit"
))) stop("Unexpected creatinine-series columns.")

if (!identical(names(episodes), c(
  "source", "episode_id", "patient_id", "t0", "age", "sex", "asa",
  "emergency", "hosp_los_days", "icu_los_days", "discharge_time", "death_time"
))) stop("Unexpected episode-covariate columns.")

if (nrow(series) && any(series$t_days < -LOOKBACK_DAYS | series$t_days > FOLLOWUP_DAYS)) {
  stop("A series row is outside the requested time window.")
}
if (nrow(series) && !all(series$unit == "mg/dL")) stop("Output units are inconsistent.")
if (nrow(series) && !all(series$episode_id %chin% episodes$episode_id)) {
  stop("A series episode_id is absent from the covariate file.")
}

n_episode_with_lab <- uniqueN(series$episode_id)
n_episode_without_lab <- nrow(episodes) - n_episode_with_lab
n_cross_admission <- ser[lab_case_id != episode_id, .N]

cat(sprintf(
  "  episodes with >=1 result: %s of %s (%.1f%%)\n",
  format(n_episode_with_lab, big.mark = ","),
  format(nrow(episodes), big.mark = ","),
  100 * n_episode_with_lab / nrow(episodes)
))
cat(sprintf(
  "  episodes with no result (retained in covariates): %s\n",
  format(n_episode_without_lab, big.mark = ",")
))
cat(sprintf(
  "  cross-admission result-to-episode links: %s\n",
  format(n_cross_admission, big.mark = ",")
))
cat(sprintf(
  "  results converted from micromol/L: %s\n",
  format(ser[converted_from_umol == TRUE, .N], big.mark = ",")
))

if (nrow(series)) {
  cat(sprintf(
    "  creatinine value, mg/dL: median %.2f (IQR %.2f-%.2f), range %.2f-%.2f\n",
    median(series$value),
    quantile(series$value, 0.25),
    quantile(series$value, 0.75),
    min(series$value),
    max(series$value)
  ))
  cat(sprintf(
    "  timing: pre-t0 %.1f%% | 0-7d %.1f%% | >7-90d %.1f%%\n",
    100 * mean(series$t_days < 0),
    100 * mean(series$t_days >= 0 & series$t_days <= 7),
    100 * mean(series$t_days > 7 & series$t_days <= 90)
  ))
}

cat(sprintf("  output: %s\n", normalizePath(series_path, winslash = "/")))
cat(sprintf("  output: %s\n", normalizePath(covariate_path, winslash = "/")))
cat("\nDone. Send both CSV files together.\n")
