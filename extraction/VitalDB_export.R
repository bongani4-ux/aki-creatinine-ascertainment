# =============================================================================
# VitalDB  -  creatinine series + case-level covariates
#
# FILES NEEDED from the PhysioNet VitalDB 1.0.0 download (only these two):
#     clinical_data.csv   one row per case (operation)
#     lab_data.csv        caseid, dt, name, result
#
# Produces (in out_dir):
#     VitalDB_creatinine_series.csv
#     VitalDB_episode_covariates.csv
#
# TIME AXIS: every time in VitalDB is seconds from that case's casestart, and
# lab_data is keyed by caseid, so a result can only be placed on the time axis
# of its own case: t_days = (dt - opstart) / 86400. lab_data covers about 90
# days either side of the case, so baseline tiers 3-4 (91-365 d) cannot occur.
#
# CREATININE: lab_data name "cr" (lab_parameters.csv: Creatinine, mg/dL,
# reference 0.70-1.40). "ccr" is creatinine clearance and is excluded.
# VitalDB records no death timestamp; death_inhosp is kept as a flag only.
# =============================================================================

rm(list = ls())
library(data.table)

# ---------------------------------------------------------------- parameters
vitaldb_dir <- "path/to/vitaldb-a-high-fidelity-multi-parameter-vital-signs-database-in-surgical-patients-1.0.0"
out_dir     <- "path/to/export_series/"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90
SEC_PER_DAY   <- 86400
CREAT_NAME    <- "cr"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
find_tbl <- function(dir, file) {
  hits <- list.files(dir, pattern = paste0("^", file, "$"), recursive = TRUE, full.names = TRUE)
  hits <- hits[!dir.exists(hits)]
  if (!length(hits)) stop("Could not find ", file, " under: ", dir)
  cat(sprintf("  %-18s -> %s\n", file, hits[1])); hits[1]
}
f_clin <- find_tbl(vitaldb_dir, "clinical_data.csv")
f_lab  <- find_tbl(vitaldb_dir, "lab_data.csv")

num <- function(x) suppressWarnings(as.numeric(x))

# ---------------------------------------------------------------- 1. episodes
clin <- fread(f_clin)
need <- c("caseid", "subjectid", "opstart", "dis", "age", "sex", "asa", "emop",
          "icu_days", "death_inhosp")
miss <- setdiff(need, names(clin))
if (length(miss)) stop("clinical_data.csv missing: ", paste(miss, collapse = ", "))

episodes <- clin[, .(
  source         = "VitalDB",
  episode_id     = as.character(caseid),
  patient_id     = as.character(subjectid),
  t0             = num(opstart),              # seconds from casestart
  age            = num(age),
  sex            = as.character(sex),
  asa            = num(asa),
  emergency      = num(emop),
  icu_los_days   = num(icu_days),
  discharge_time = num(dis),                  # seconds from casestart
  discharge_t_days = (num(dis) - num(opstart)) / SEC_PER_DAY,
  death_time     = NA_real_,                  # not recorded in VitalDB
  died_inhosp    = num(death_inhosp)
)]
fwrite(episodes, file.path(out_dir, "VitalDB_episode_covariates.csv"))
cat(sprintf("  cases written: %s (%s subjects)\n", format(nrow(episodes), big.mark = ","),
            format(uniqueN(episodes$patient_id), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
lab <- fread(f_lab)
cat("  lab names containing 'cr':", paste(sort(unique(grep("cr", lab$name, value = TRUE))),
                                          collapse = ", "), "\n")
cr <- lab[name == CREAT_NAME, .(episode_id = as.character(caseid), dt = num(dt),
                                value = num(result))]
cr <- cr[!is.na(dt) & !is.na(value) & value > 0 & value < 50]
if (!nrow(cr)) stop("No rows with name == '", CREAT_NAME, "'. Check the printed lab names.")

ser <- merge(cr, episodes[, .(episode_id, patient_id, t0)], by = "episode_id")
ser[, t_days := (dt - t0) / SEC_PER_DAY]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]
series <- ser[, .(source = "VitalDB", episode_id, patient_id, timestamp = dt,
                  t_days = round(t_days, 5), value, unit = "mg/dL")]
setorder(series, episode_id, t_days)
fwrite(series, file.path(out_dir, "VitalDB_creatinine_series.csv"))

# ---------------------------------------------------------------- 3. checks
m <- median(series$value)
cat(sprintf("  series written: %s rows | median %.2f mg/dL (expect ~0.7-1.0)\n",
            format(nrow(series), big.mark = ","), m))
if (m > 20) stop("Values look like umol/L, not mg/dL; check lab_parameters.csv.")
cat(sprintf("  cases with >=1 result: %.1f%% | pre-t0 results %.1f%%\n",
            100 * uniqueN(series$episode_id) / nrow(episodes), 100 * mean(series$t_days < 0)))
cat("Done. Do not open the output in Excel.\n")
