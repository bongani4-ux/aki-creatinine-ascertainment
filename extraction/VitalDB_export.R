# clinical_data and lab_data are already loaded; do not clear the workspace.
library(data.table)

# Assumes clinical_data and lab_data are already loaded in the R session.
# If dt is in seconds, leave TIME_UNIT = "secs"; change to "mins" if needed.
LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90
TIME_UNIT <- "secs"
LAB_UNIT <- "umol/L"       # change to "mg/dL" if lab_data$result is already mg/dL
OUT_DIR <- "path/to/vitaldb"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

stopifnot(exists("clinical_data"), exists("lab_data"))
clinical <- as.data.table(copy(clinical_data))
labs <- as.data.table(copy(lab_data))

required_clin <- c("caseid", "subjectid", "opstart", "dis", "icu_days",
                   "age", "sex", "asa", "emop", "death_inhosp")
required_lab <- c("caseid", "dt", "name", "result")
miss_clin <- setdiff(required_clin, names(clinical))
miss_lab <- setdiff(required_lab, names(labs))
if (length(miss_clin)) stop("clinical_data missing: ", paste(miss_clin, collapse = ", "))
if (length(miss_lab)) stop("lab_data missing: ", paste(miss_lab, collapse = ", "))

num <- function(x) suppressWarnings(as.numeric(x))
id <- function(x) as.character(x)

# VitalDB has no death timestamp in the supplied clinical schema.
# death_time is therefore NA; death_inhosp is retained as a flag.
episodes <- clinical[, .(
  episode_id = id(caseid),
  patient_id = id(subjectid),
  t0 = num(opstart),
  age = num(age),
  sex = as.character(sex),
  ASA = num(asa),
  emergency = num(emop),
  hosp_los_days = num(dis) / if (TIME_UNIT == "secs") 86400 else 1440,
  icu_los_days = num(icu_days),
  discharge_time = num(dis),
  death_time = NA_real_,
  death_inhosp = num(death_inhosp)
)]

fwrite(episodes, file.path(OUT_DIR, "episode_covariates.csv"))

# Keep every lab row; filtering by name is case-insensitive and accepts common aliases.
cr <- labs[tolower(trimws(as.character(name))) %chin%
             c("creatinine", "serum creatinine", "creat") &
             !is.na(result) & num(result) > 0,
           .(episode_id = id(caseid), timestamp_raw = num(dt), value_raw = num(result))]
cr <- cr[!is.na(timestamp_raw) & !is.na(value_raw)]

# The supplied lab_data is keyed by caseid, so this join preserves available
# indexed-case labs. If a patient-level lab table becomes available, replace
# episode_id with patient_id and join every result to every episode of that patient.
series <- merge(cr, episodes[, .(episode_id, patient_id, t0)],
                 by = "episode_id", all.x = TRUE)
series[, t_days := (timestamp_raw - t0) /
         if (TIME_UNIT == "secs") 86400 else 1440]
series <- series[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

if (LAB_UNIT == "umol/L") {
  series[, value := value_raw / 88.4]
  output_unit <- "mg/dL"
  cat("Converted creatinine from umol/L to mg/dL using value / 88.4.\n")
} else {
  series[, value := value_raw]
  output_unit <- "mg/dL"
  cat("Creatinine treated as already in mg/dL; no conversion applied.\n")
}

series <- series[, .(episode_id, patient_id,
                     timestamp = timestamp_raw,
                     t_days = round(t_days, 5), value,
                     unit = output_unit)]
setorder(series, episode_id, t_days)
fwrite(series, file.path(OUT_DIR, "creatinine_series.csv"))

cat(sprintf("Wrote %s episodes and %s creatinine results.\n",
            format(nrow(episodes), big.mark = ","),
            format(nrow(series), big.mark = ",")))
cat("Use data.table::fread() to inspect the CSV files; do not open them in Excel.\n")
