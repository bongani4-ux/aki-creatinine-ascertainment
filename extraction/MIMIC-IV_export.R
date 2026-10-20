# =============================================================================
# MIMIC-IV
# Full creatinine measurement series + episode-level covariates
#
# Produces two files:
#   MIMIC-IV_creatinine_series.csv      one row per creatinine result
#   MIMIC-IV_episode_covariates.csv     one row per ICU stay
#
# Satisfies reporting-set item 7 (distribution of number and timing of
# creatinine measurements per patient) and supplies age, sex, LOS and mortality.
# =============================================================================

rm(list = ls())
library(data.table)
library(duckdb)
library(DBI)

# ---------------------------------------------------------------- parameters
out_dir <- "path/to/export_series/"

creat_itemid  <- 50912    # serum creatinine, already mg/dL in MIMIC-IV
LOOKBACK_DAYS <- 365      # export window before t0
FOLLOWUP_DAYS <- 90       # export window after t0

# If you already unzipped for the earlier corrected extraction, reuse that folder.
extract_dir <- "path/to/mimic_extracted"
zip_path    <- "path/to/mimic-iv-3.1.zip"

get_file <- function(internal, dest_dir, zip) {
  f <- file.path(dest_dir, basename(internal))
  if (file.exists(f)) { cat("  have:", basename(f), "\n"); return(f) }
  cat("  extracting:", basename(internal), "\n")
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  unzip(zip, files = internal, exdir = dest_dir, junkpaths = TRUE)
  if (!file.exists(f)) stop("Not found in zip: ", internal)
  f
}

cat("=== Locating files ===\n")
f_icustays   <- get_file("mimic-iv-3.1/icu/icustays.csv.gz",   extract_dir, zip_path)
f_patients   <- get_file("mimic-iv-3.1/hosp/patients.csv.gz",  extract_dir, zip_path)
f_admissions <- get_file("mimic-iv-3.1/hosp/admissions.csv.gz",extract_dir, zip_path)
f_labevents  <- get_file("mimic-iv-3.1/hosp/labevents.csv.gz", extract_dir, zip_path)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------------------------------------------- 1. episodes
cat("\n=== Building episode table ===\n")
icu <- fread(f_icustays)
icu[, `:=`(stay_id = as.character(stay_id),
           subject_id = as.character(subject_id),
           hadm_id = as.character(hadm_id),
           intime  = as.POSIXct(intime,  tz = "UTC"),
           outtime = as.POSIXct(outtime, tz = "UTC"))]

pat <- fread(f_patients)[, .(subject_id = as.character(subject_id),
                             gender, anchor_age, dod)]
adm <- fread(f_admissions)[, .(hadm_id = as.character(hadm_id),
                               admittime, dischtime, deathtime,
                               hospital_expire_flag, admission_type, race)]
adm[, `:=`(admittime = as.POSIXct(admittime, tz = "UTC"),
           dischtime = as.POSIXct(dischtime, tz = "UTC"),
           deathtime = as.POSIXct(deathtime, tz = "UTC"))]

ep <- merge(icu, pat, by = "subject_id", all.x = TRUE)
ep <- merge(ep, adm, by = "hadm_id", all.x = TRUE)

ep[, `:=`(
  t0            = intime,                                   # time zero (C1)
  age           = anchor_age,
  sex           = gender,
  icu_los_days  = as.numeric(difftime(outtime, intime, units = "days")),
  hosp_los_days = as.numeric(difftime(dischtime, admittime, units = "days")),
  died_inhosp   = as.integer(hospital_expire_flag == 1)
)]
# days from t0 to death (NA if survived); Inf-safe
ep[, death_t_days := as.numeric(difftime(deathtime, intime, units = "days"))]
ep[is.na(death_t_days) & !is.na(dod),
   death_t_days := as.numeric(difftime(as.POSIXct(dod, tz = "UTC"), intime, units = "days"))]

episodes <- ep[, .(
  source = "MIMIC-IV",
  episode_id = stay_id,
  patient_id = subject_id,
  hadm_id,
  t0,
  age, sex, race,
  admission_type,
  icu_los_days, hosp_los_days,
  died_inhosp,
  death_t_days,
  discharge_t_days = as.numeric(difftime(dischtime, intime, units = "days"))
)]

fwrite(episodes, paste0(out_dir, "MIMIC-IV_episode_covariates.csv"))
cat(sprintf("  episodes: %s rows written\n", format(nrow(episodes), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
cat("\n=== Extracting creatinine series (this takes several minutes) ===\n")
con <- dbConnect(duckdb(), dbdir = ":memory:")
q <- sprintf("
  SELECT subject_id, hadm_id, charttime, valuenum
  FROM read_csv_auto('%s', header=true, strict_mode=false,
                     ignore_errors=true, quote='\"', escape='\"')
  WHERE itemid = %d AND valuenum IS NOT NULL
    AND valuenum > 0 AND valuenum < 50
", gsub("\\\\", "/", f_labevents), creat_itemid)
cr <- as.data.table(dbGetQuery(con, q))
dbDisconnect(con, shutdown = TRUE)
cat(sprintf("  raw creatinine rows: %s\n", format(nrow(cr), big.mark = ",")))

cr[, `:=`(subject_id = as.character(subject_id),
          charttime  = as.POSIXct(charttime, tz = "UTC"))]

# Join to EVERY ICU stay of that subject, then window on t_days.
# A subject with 2 stays contributes each result twice, once per episode —
# this is correct: each episode has its own t0.
ser <- merge(cr[, .(subject_id, charttime, value = valuenum)],
             episodes[, .(episode_id, patient_id, t0)],
             by.x = "subject_id", by.y = "patient_id",
             allow.cartesian = TRUE)
ser[, t_days := as.numeric(difftime(charttime, t0, units = "days"))]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

series <- ser[, .(
  source     = "MIMIC-IV",
  episode_id,
  patient_id = subject_id,
  charttime,
  t_days     = round(t_days, 5),
  value,
  unit       = "mg/dL"
)]
setorder(series, episode_id, t_days)

fwrite(series, paste0(out_dir, "MIMIC-IV_creatinine_series.csv"))
cat(sprintf("  series: %s rows written\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
cat("\n=== Sanity checks ===\n")
cat(sprintf("  median creatinine: %.2f mg/dL  (expect ~0.8-1.1)\n", median(series$value)))
cat(sprintf("  episodes with >=1 result: %s of %s (%.1f%%)\n",
            format(uniqueN(series$episode_id), big.mark = ","),
            format(nrow(episodes), big.mark = ","),
            100 * uniqueN(series$episode_id) / nrow(episodes)))
n_per <- series[, .N, by = episode_id]
cat(sprintf("  results per episode: median %d (IQR %d-%d), max %d\n",
            as.integer(median(n_per$N)),
            as.integer(quantile(n_per$N, .25)),
            as.integer(quantile(n_per$N, .75)),
            max(n_per$N)))
cat(sprintf("  pre-t0 results:  %s (%.1f%%)\n",
            format(sum(series$t_days < 0), big.mark = ","),
            100 * mean(series$t_days < 0)))
cat(sprintf("  day 7-90 results: %s (%.1f%%)\n",
            format(sum(series$t_days > 7 & series$t_days <= 90), big.mark = ","),
            100 * mean(series$t_days > 7 & series$t_days <= 90)))
cat(sprintf("  age non-missing: %.1f%%   sex non-missing: %.1f%%\n",
            100 * mean(!is.na(episodes$age)), 100 * mean(!is.na(episodes$sex))))

cat("\nDone. Upload BOTH files:\n")
cat("  ", paste0(out_dir, "MIMIC-IV_creatinine_series.csv"), "\n")
cat("  ", paste0(out_dir, "MIMIC-IV_episode_covariates.csv"), "\n")
cat("If the series file is very large, split by episode_id or gzip it before upload.\n")
