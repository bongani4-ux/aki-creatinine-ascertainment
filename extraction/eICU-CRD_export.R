suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table first: install.packages('data.table')", call. = FALSE)
  }
  library(data.table)
})

# Export the eICU creatinine time series and one-row-per-ICU-stay covariates.
#
# Usage:
#   Rscript eICU-CRD_export.R path/to/eicu-collaborative-research-database-2.0 path/to/export_series
# Writes eICU-CRD_creatinine_series.csv and eICU-CRD_episode_covariates.csv.
#
# eICU stores event times as minute offsets rather than full dates. Accordingly:
#   * episode_id is patientunitstayid and patient_id is uniquepid.
#   * t0 is ICU admission and is represented as 0 minutes.
#   * timestamp, discharge_time, and death_time are minute offsets from t0.
#   * t_days is timestamp / 1440.
#
# Lab rows from every ICU unit stay in the same patienthealthsystemstayid are
# aligned to each index ICU stay. Separate hospitalizations and outpatient data
# cannot be aligned because eICU does not provide the required cross-encounter
# dates; the requested -365 to +90 day filter is still applied to all alignable
# rows.
#
# First episode per patient: every stay has t0 = 0, so the pipeline cannot order a
# patient's stays in time and falls back to episode_id order (docs/methods_notes.md).
# unitvisitnumber and hospitaldischargeyear are exported to allow a time-ordered
# sensitivity analysis.

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else "path/to/eicu-collaborative-research-database-2.0"
out_dir <- if (length(args) >= 2L) args[[2L]] else "path/to/export_series"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90
MINUTES_PER_DAY <- 1440
UMOL_PER_MG_DL <- 88.4

patient_path <- file.path(data_dir, "patient.csv")
lab_path <- file.path(data_dir, "lab.csv")
if (!file.exists(patient_path)) stop("Missing input: ", patient_path)
if (!file.exists(lab_path)) stop("Missing input: ", lab_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_names <- function(x) sub("^\ufeff", "", x)
blank_or_na <- function(x) is.na(x) | trimws(as.character(x)) == ""
require_columns <- function(actual, expected, table_name) {
  missing <- setdiff(expected, actual)
  if (length(missing)) {
    stop(table_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

message("Reading patient.csv (only required fields) ...")
patient_cols <- c(
  "patientunitstayid", "patienthealthsystemstayid", "uniquepid", "gender", "age",
  "hospitaladmitoffset", "hospitaladmitsource", "hospitaldischargeoffset",
  "hospitaldischargestatus", "unitdischargeoffset", "unitvisitnumber",
  "hospitaldischargeyear"
)
patients <- fread(patient_path, select = patient_cols, showProgress = TRUE)
setnames(patients, clean_names(names(patients)))
require_columns(names(patients), patient_cols, "patient.csv")

patients[, `:=`(
  patientunitstayid = as.character(patientunitstayid),
  patienthealthsystemstayid = as.character(patienthealthsystemstayid),
  uniquepid = as.character(uniquepid),
  hospitaladmitoffset = suppressWarnings(as.numeric(hospitaladmitoffset)),
  hospitaldischargeoffset = suppressWarnings(as.numeric(hospitaldischargeoffset)),
  unitdischargeoffset = suppressWarnings(as.numeric(unitdischargeoffset))
)]

if (anyDuplicated(patients$patientunitstayid)) {
  stop("patientunitstayid is not unique in patient.csv; episode grain is invalid.")
}
if (any(blank_or_na(patients$patientunitstayid)) ||
    any(blank_or_na(patients$patienthealthsystemstayid)) ||
    any(blank_or_na(patients$uniquepid))) {
  stop("Required eICU identifiers contain missing values.")
}

patient_conflicts <- patients[, .(n_patient_ids = uniqueN(uniquepid)),
                              by = patienthealthsystemstayid][n_patient_ids != 1L]
if (nrow(patient_conflicts)) {
  stop("A patienthealthsystemstayid maps to more than one uniquepid in patient.csv.")
}

# eICU masks age above 89 as the string "> 89". Use 90 as the documented
# numeric surrogate rather than treating this protected group as missing.
age_text <- trimws(as.character(patients$age))
age_num <- suppressWarnings(as.numeric(age_text))
age_num[age_text == "> 89"] <- 90
age_num[age_num < 0 | age_num > 110] <- NA_real_

sex_text <- tolower(trimws(as.character(patients$gender)))
sex <- fcase(
  sex_text == "male", "M",
  sex_text == "female", "F",
  default = NA_character_
)

# eICU has admission source rather than an elective/emergency flag. This is a
# reproducible proxy: 1 only for hospital admission through the ED, 0 for another
# recorded source, and missing when hospital admission source is blank.
admit_source <- trimws(as.character(patients$hospitaladmitsource))
emergency <- fcase(
  admit_source == "Emergency Department", 1L,
  !blank_or_na(admit_source), 0L,
  default = NA_integer_
)

hosp_los_days <- (patients$hospitaldischargeoffset - patients$hospitaladmitoffset) /
  MINUTES_PER_DAY
icu_los_days <- patients$unitdischargeoffset / MINUTES_PER_DAY
discharge_time <- patients$hospitaldischargeoffset
death_time <- fifelse(
  trimws(as.character(patients$hospitaldischargestatus)) == "Expired",
  patients$hospitaldischargeoffset,
  NA_real_
)

# Negative lengths of stay and discharge/death offsets are temporally impossible
# for an ICU-admission anchor. Preserve the episodes but represent these small
# source-data inconsistencies as missing covariates.
invalid_hosp_los_rows <- sum(!is.na(hosp_los_days) & hosp_los_days < 0)
invalid_icu_los_rows <- sum(!is.na(icu_los_days) & icu_los_days < 0)
invalid_discharge_time_rows <- sum(!is.na(discharge_time) & discharge_time < 0)
invalid_death_time_rows <- sum(!is.na(death_time) & death_time < 0)
hosp_los_days[!is.na(hosp_los_days) & hosp_los_days < 0] <- NA_real_
icu_los_days[!is.na(icu_los_days) & icu_los_days < 0] <- NA_real_
discharge_time[!is.na(discharge_time) & discharge_time < 0] <- NA_real_
death_time[!is.na(death_time) & death_time < 0] <- NA_real_

episode_covariates <- data.table(
  episode_id = patients$patientunitstayid,
  patient_id = patients$uniquepid,
  t0 = 0,
  age = age_num,
  sex = sex,
  ASA = NA_real_,
  emergency = emergency,
  hosp_los_days = hosp_los_days,
  icu_los_days = icu_los_days,
  discharge_time = discharge_time,
  death_time = death_time,
  unitvisitnumber = patients$unitvisitnumber,
  hospitaldischargeyear = patients$hospitaldischargeyear
)

expected_covariate_columns <- c(
  "episode_id", "patient_id", "t0", "age", "sex", "ASA", "emergency",
  "hosp_los_days", "icu_los_days", "discharge_time", "death_time",
  "unitvisitnumber", "hospitaldischargeyear"
)
stopifnot(identical(names(episode_covariates), expected_covariate_columns))
stopifnot(nrow(episode_covariates) == nrow(patients))
stopifnot(!anyDuplicated(episode_covariates$episode_id))

covariate_path <- file.path(out_dir, "eICU-CRD_episode_covariates.csv")
fwrite(episode_covariates, covariate_path, na = "")
message("Wrote ", format(nrow(episode_covariates), big.mark = ","),
        " episode rows to ", covariate_path)

message("Reading lab.csv (only creatinine-relevant fields) ...")
lab_cols <- c(
  "patientunitstayid", "labresultoffset", "labname", "labresult",
  "labmeasurenamesystem", "labmeasurenameinterface"
)
labs <- fread(lab_path, select = lab_cols, showProgress = TRUE)
setnames(labs, clean_names(names(labs)))
require_columns(names(labs), lab_cols, "lab.csv")

creatinine_label <- tolower(trimws(as.character(labs$labname))) == "creatinine"
cr <- labs[creatinine_label, .(
  lab_stay_id = as.character(patientunitstayid),
  labresultoffset = suppressWarnings(as.numeric(labresultoffset)),
  value_raw = suppressWarnings(as.numeric(labresult)),
  unit_system = trimws(as.character(labmeasurenamesystem)),
  unit_interface = trimws(as.character(labmeasurenameinterface))
)]
rm(labs)
invisible(gc())

if (!nrow(cr)) stop("No lab rows with labname == 'creatinine' were found.")

cr[, unit_raw := fifelse(!blank_or_na(unit_system), unit_system, unit_interface)]
cr[, unit_norm := tolower(gsub("\\s+", "", unit_raw))]

# eICU normally stores creatinine in mg/dL. The conversion branch makes the
# export safe if a local copy contains micromolar values.
is_mg_dl <- blank_or_na(cr$unit_norm) | cr$unit_norm %chin% c(
  "mg/dl", "mgdl", "mg/dl."
)
is_umol_l <- grepl("umol/l|umolperliter|micromol/l|micromole/l|µmol/l|μmol/l",
                   cr$unit_norm, perl = TRUE)
unsupported_units <- unique(cr$unit_raw[!is_mg_dl & !is_umol_l])
unsupported_units <- unsupported_units[!blank_or_na(unsupported_units)]
if (length(unsupported_units)) {
  stop("Unsupported creatinine unit(s): ", paste(unsupported_units, collapse = ", "))
}

cr[, value := fifelse(is_umol_l, value_raw / UMOL_PER_MG_DL, value_raw)]
converted_rows <- sum(is_umol_l, na.rm = TRUE)
invalid_value_rows <- cr[is.na(labresultoffset) | is.na(value) | value <= 0 | value >= 50, .N]
cr <- cr[!is.na(labresultoffset) & !is.na(value) & value > 0 & value < 50]
if (!nrow(cr)) stop("No valid creatinine rows remained after value/time validation.")

lab_stay_map <- patients[, .(
  lab_stay_id = patientunitstayid,
  patienthealthsystemstayid,
  lab_hospitaladmitoffset = hospitaladmitoffset
)]
cr <- lab_stay_map[cr, on = "lab_stay_id", nomatch = 0L]
cr[, lab_from_hospital_admit_min := labresultoffset - lab_hospitaladmitoffset]

episode_map <- patients[, .(
  patienthealthsystemstayid,
  episode_id = patientunitstayid,
  patient_id = uniquepid,
  episode_hospitaladmitoffset = hospitaladmitoffset
)]

# A lab event is paired to every ICU episode within the same hospitalization.
# Its episode-relative timestamp equals the lab's hospital-relative offset minus
# the index ICU admission's hospital-relative offset.
series_work <- merge(
  cr,
  episode_map,
  by = "patienthealthsystemstayid",
  allow.cartesian = TRUE,
  sort = FALSE
)
series_work[, timestamp := lab_from_hospital_admit_min + episode_hospitaladmitoffset]
series_work[, t_days := timestamp / MINUTES_PER_DAY]
series_work <- series_work[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

cross_unit_rows <- series_work[lab_stay_id != episode_id, .N]
creatinine_series <- series_work[, .(
  episode_id,
  patient_id,
  timestamp,
  t_days = round(t_days, 5),
  value = round(value, 6),
  unit = "mg/dL"
)]
setorder(creatinine_series, episode_id, timestamp)

expected_series_columns <- c(
  "episode_id", "patient_id", "timestamp", "t_days", "value", "unit"
)
stopifnot(identical(names(creatinine_series), expected_series_columns))
stopifnot(all(creatinine_series$episode_id %chin% episode_covariates$episode_id))
stopifnot(all(creatinine_series$t_days >= -LOOKBACK_DAYS &
              creatinine_series$t_days <= FOLLOWUP_DAYS))
stopifnot(all(creatinine_series$unit == "mg/dL"))
stopifnot(max(abs(creatinine_series$t_days -
                  creatinine_series$timestamp / MINUTES_PER_DAY)) <= 0.00001)

series_path <- file.path(out_dir, "eICU-CRD_creatinine_series.csv")
fwrite(creatinine_series, series_path, na = "")

n_with_creatinine <- uniqueN(creatinine_series$episode_id)
n_without_creatinine <- nrow(episode_covariates) - n_with_creatinine
n_exact_duplicates <- sum(duplicated(creatinine_series))
message("\n=== Export quality summary ===")
message("Episode grain: one row per patientunitstayid")
message("Episodes: ", format(nrow(episode_covariates), big.mark = ","))
message("Unique patients: ", format(uniqueN(episode_covariates$patient_id), big.mark = ","))
message("Creatinine series rows: ", format(nrow(creatinine_series), big.mark = ","))
message("Episodes with >=1 creatinine: ", format(n_with_creatinine, big.mark = ","),
        " (", sprintf("%.2f%%", 100 * n_with_creatinine / nrow(episode_covariates)), ")")
message("Episodes with no creatinine retained in covariates: ",
        format(n_without_creatinine, big.mark = ","))
message("Rows contributed by another ICU unit in the same hospitalization: ",
        format(cross_unit_rows, big.mark = ","))
message("Creatinine rows converted from umol/L to mg/dL: ",
        format(converted_rows, big.mark = ","))
message("Invalid raw creatinine rows excluded: ",
        format(invalid_value_rows, big.mark = ","))
message("Median creatinine: ", sprintf("%.3f mg/dL", median(creatinine_series$value)))
message("Exact duplicate output rows retained: ", format(n_exact_duplicates, big.mark = ","))
message("Invalid negative hospital LOS values set to missing: ",
        format(invalid_hosp_los_rows, big.mark = ","))
message("Invalid negative ICU LOS values set to missing: ",
        format(invalid_icu_los_rows, big.mark = ","))
message("Invalid negative discharge offsets set to missing: ",
        format(invalid_discharge_time_rows, big.mark = ","))
message("Invalid negative death offsets set to missing: ",
        format(invalid_death_time_rows, big.mark = ","))
message("Age non-missing: ", sprintf("%.2f%%", 100 * mean(!is.na(episode_covariates$age))))
message("Sex non-missing: ", sprintf("%.2f%%", 100 * mean(!is.na(episode_covariates$sex))))
message("Emergency proxy non-missing: ",
        sprintf("%.2f%%", 100 * mean(!is.na(episode_covariates$emergency))))
message("ASA non-missing: 0.00% (field is unavailable in eICU)")
message("Time range: ", sprintf("%.5f", min(creatinine_series$t_days)),
        " to ", sprintf("%.5f", max(creatinine_series$t_days)), " days")
message("\nWrote: ", series_path)
message("Wrote: ", covariate_path)
message("No Excel files were created.")
