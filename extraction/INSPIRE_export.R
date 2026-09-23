# =============================================================================
# EXPORT 2 of 3 — INSPIRE  (v1.4.2 compatible)
# Full creatinine measurement series + operation-level covariates
#
# FILES NEEDED from the PhysioNet INSPIRE download — only two:
#     operations.csv (or operations.csv.gz)
#     labs.csv       (or labs.csv.gz)
# You do NOT need vitals, ward_vitals, diagnosis, medications or parameters.
#
# Produces:
#     INSPIRE_creatinine_series.csv
#     INSPIRE_operation_covariates.csv
#
# TIME AXIS: INSPIRE timestamps are integers (minutes) offset from a
# de-identified per-subject origin. opdate and lab chart_time share that axis,
# so t_days = (chart_time - opdate) / 1440.
# =============================================================================

rm(list = ls())
library(data.table)

# ---------------------------------------------------------------- parameters
inspire_dir <- "path/to/inspire-a-publicly-available-research-dataset-for-perioperative-medicine-1.4.2/inspire-a-publicly-available-research-dataset-for-perioperative-medicine-1.4.2"
out_dir     <- "path/to/export_series/"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90
MIN_PER_DAY   <- 1440

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------------------------------------------- path finder
# In the 1.4.2 release each table appears more than once at the top level:
#   operations.csv     -> a FOLDER that itself contains operations.csv
#   operations.csv.gz  -> a real gzipped file
# This resolves whichever exists, preferring the .gz file.
resolve <- function(dir, stem) {
  cands <- c(
    file.path(dir, paste0(stem, ".csv.gz")),
    file.path(dir, stem, paste0(stem, ".csv")),
    file.path(dir, paste0(stem, ".csv"))
  )
  for (p in cands) if (file.exists(p) && !dir.exists(p)) {
    cat(sprintf("  %-12s -> %s\n", stem, p)); return(p)
  }
  stop("Could not find a readable file for '", stem, "' under:\n  ", dir,
       "\nTried:\n  ", paste(cands, collapse = "\n  "))
}

cat("=== Resolving input files ===\n")
f_ops  <- resolve(inspire_dir, "operations")
f_labs <- resolve(inspire_dir, "labs")

# ---------------------------------------------------------------- 0. inspect
cat("\n=== Schema check (read this before trusting the output) ===\n")
ops_head <- fread(f_ops, nrows = 5)
cat("operations columns:\n"); print(names(ops_head))

lab_head <- fread(f_labs, nrows = 200000)
setnames(lab_head, tolower(names(lab_head)))
cat("\nlabs columns:\n"); print(names(lab_head))

item_col_h <- intersect(c("item_name","item","label","test_name","itemid"), names(lab_head))[1]
cat(sprintf("\nDistinct lab labels containing 'creat' (column '%s'):\n", item_col_h))
creat_labels <- unique(grep("creat", lab_head[[item_col_h]], ignore.case = TRUE, value = TRUE))
print(creat_labels)
if (length(creat_labels) == 0) {
  cat("\n!! No label matched 'creat'. First 40 distinct labels in this file:\n")
  print(head(unique(lab_head[[item_col_h]]), 40))
  stop("Edit CREAT_LABEL below to match the correct label, then re-run.")
}

CREAT_LABEL <- "creatinine"   # <-- edit if the printed list shows a different label

# ---------------------------------------------------------------- 1. episodes
cat("\n=== Building operation table ===\n")
ops <- fread(f_ops)
setnames(ops, tolower(names(ops)))

g <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else rep(NA, nrow(dt))
mins_to_days <- function(x, t0) (suppressWarnings(as.numeric(x)) - t0) / MIN_PER_DAY

t0 <- suppressWarnings(as.numeric(ops$opdate))

episodes <- data.table(
  source       = "INSPIRE",
  episode_id   = as.character(ops$op_id),
  patient_id   = as.character(ops$subject_id),
  hadm_id      = as.character(g(ops, "hadm_id")),
  t0_min       = t0,
  age          = suppressWarnings(as.numeric(g(ops, "age"))),
  sex          = as.character(g(ops, "sex")),
  weight       = suppressWarnings(as.numeric(g(ops, "weight"))),
  height       = suppressWarnings(as.numeric(g(ops, "height"))),
  asa          = suppressWarnings(as.numeric(g(ops, "asa"))),
  emergency    = suppressWarnings(as.numeric(g(ops, "emop"))),
  department   = as.character(g(ops, "department")),
  anaes_type   = as.character(g(ops, "antype")),
  icd10_pcs    = as.character(g(ops, "icd10_pcs"))
)

episodes[, `:=`(
  admission_t_days      = mins_to_days(g(ops, "admission_time"),      t0),
  discharge_t_days      = mins_to_days(g(ops, "discharge_time"),      t0),
  opstart_t_days        = mins_to_days(g(ops, "opstart_time"),        t0),
  opend_t_days          = mins_to_days(g(ops, "opend_time"),          t0),
  icu_in_t_days         = mins_to_days(g(ops, "icuin_time"),          t0),
  icu_out_t_days        = mins_to_days(g(ops, "icuout_time"),         t0),
  death_t_days          = mins_to_days(g(ops, "inhosp_death_time"),   t0),
  allcause_death_t_days = mins_to_days(g(ops, "allcause_death_time"), t0)
)]
episodes[, `:=`(
  hosp_los_days = discharge_t_days - admission_t_days,
  icu_los_days  = icu_out_t_days - icu_in_t_days,
  icu_admitted  = as.integer(!is.na(icu_in_t_days)),
  died_inhosp   = as.integer(!is.na(death_t_days))
)]
episodes[age < 0 | age > 110, age := NA]
episodes[hosp_los_days < 0 | hosp_los_days > 365, hosp_los_days := NA]

fwrite(episodes, paste0(out_dir, "INSPIRE_operation_covariates.csv"))
cat(sprintf("  operations written: %s rows\n", format(nrow(episodes), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
cat("\n=== Extracting creatinine series ===\n")
labs <- fread(f_labs)
setnames(labs, tolower(names(labs)))

item_col  <- intersect(c("item_name","item","label","test_name"), names(labs))[1]
time_col  <- intersect(c("chart_time","charttime","time","lab_time"), names(labs))[1]
value_col <- intersect(c("value","valuenum","result","result_value"), names(labs))[1]
cat(sprintf("  using: item=%s  time=%s  value=%s\n", item_col, time_col, value_col))
if (any(is.na(c(item_col, time_col, value_col))))
  stop("Could not identify a required lab column. Names are:\n  ",
       paste(names(labs), collapse = ", "))

cr <- labs[grepl(CREAT_LABEL, get(item_col), ignore.case = TRUE),
           .(patient_id = as.character(subject_id),
             chart_min  = suppressWarnings(as.numeric(get(time_col))),
             value      = suppressWarnings(as.numeric(get(value_col))))]
cr <- cr[!is.na(value) & !is.na(chart_min) & value > 0 & value < 50]
cat(sprintf("  creatinine results found: %s\n", format(nrow(cr), big.mark = ",")))
if (nrow(cr) == 0) stop("No creatinine rows after filtering. Check CREAT_LABEL.")

# Join every result to every operation of that subject, then window on t_days.
# A subject with two operations contributes each result twice, once per episode
# — correct, because each operation has its own t0.
ser <- merge(cr, episodes[, .(episode_id, patient_id, t0_min)],
             by = "patient_id", allow.cartesian = TRUE)
ser[, t_days := (chart_min - t0_min) / MIN_PER_DAY]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

series <- ser[, .(
  source     = "INSPIRE",
  episode_id, patient_id,
  chart_min,
  t_days     = round(t_days, 5),
  value,
  unit       = "mg/dL"     # verified by the median check below
)]
setorder(series, episode_id, t_days)
fwrite(series, paste0(out_dir, "INSPIRE_creatinine_series.csv"))
cat(sprintf("  series written: %s rows\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
cat("\n=== Sanity checks ===\n")
m <- median(series$value)
cat(sprintf("  median creatinine: %.2f\n", m))
if (m > 20) {
  cat("  !! values look like umol/L, not mg/dL.\n")
  cat("     Change unit = \"mg/dL\" to \"umol/L\" in the series block and re-run.\n")
} else {
  cat("     consistent with mg/dL; unit label is correct.\n")
}
cat(sprintf("  operations with >=1 result: %s of %s (%.1f%%)\n",
            format(uniqueN(series$episode_id), big.mark = ","),
            format(nrow(episodes), big.mark = ","),
            100 * uniqueN(series$episode_id) / nrow(episodes)))
n_per <- series[, .N, by = episode_id]
cat(sprintf("  results per operation: median %d (IQR %d-%d), max %d\n",
            as.integer(median(n_per$N)),
            as.integer(quantile(n_per$N, .25)),
            as.integer(quantile(n_per$N, .75)),
            max(n_per$N)))
cat(sprintf("  timing:  pre-t0 %.1f%%  |  7-365d before %.1f%%  |  0-7d after %.1f%%  |  day 7-90 after %.1f%%\n",
            100 * mean(series$t_days < 0),
            100 * mean(series$t_days <= -7 & series$t_days >= -365),
            100 * mean(series$t_days > 0 & series$t_days <= 7),
            100 * mean(series$t_days > 7 & series$t_days <= 90)))
cat(sprintf("  covariates non-missing: age %.1f%%  sex %.1f%%  ASA %.1f%%  hosp LOS %.1f%%  death time %.1f%%\n",
            100 * mean(!is.na(episodes$age)),
            100 * mean(!is.na(episodes$sex)),
            100 * mean(!is.na(episodes$asa)),
            100 * mean(!is.na(episodes$hosp_los_days)),
            100 * mean(!is.na(episodes$death_t_days))))

cat("\nDone. Upload BOTH files:\n")
cat("  ", paste0(out_dir, "INSPIRE_creatinine_series.csv"), "\n")
cat("  ", paste0(out_dir, "INSPIRE_operation_covariates.csv"), "\n")
cat("If the series file exceeds ~200 MB, re-run the two fwrite calls with\n")
cat("compress = \"gzip\" and upload the .gz versions instead.\n")
