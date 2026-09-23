# =============================================================================
# EXPORT — PIC  (Paediatric Intensive Care database, Zhejiang University)
# Full creatinine measurement series + episode-level covariates
#
# IMPORTANT — TWO VERSIONS SHIP IN THE SAME DOWNLOAD:
#   root level   *.csv.gz          = v1.0.0
#   V1.1.0/      *.csv             = v1.1.0  <-- USE THIS ONE
# An earlier export appears to have produced only 7,488 ICU stays against the
# ~13,400 in the full dataset. This script reads V1.1.0 explicitly and prints
# the row count so the denominator can be confirmed.
#
# FILES NEEDED (all inside V1.1.0/):
#     ICUSTAYS.csv   PATIENTS.csv   ADMISSIONS.csv   LABEVENTS.csv   D_LABITEMS.csv
#
# Produces:
#     PIC_creatinine_series.csv
#     PIC_episode_covariates.csv
#
# t0 = ICU admission (INTIME). PIC is paediatric: report it separately and
# never pool it with adult sources.
# =============================================================================

rm(list = ls())
library(data.table)
suppressWarnings(suppressMessages({ library(duckdb); library(DBI) }))

# ---------------------------------------------------------------- parameters
pic_dir <- "path/to/paediatric-intensive-care-database-1.1.0/paediatric-intensive-care-database-1.1.0"
out_dir <- "path/to/export_series/"

USE_V110      <- TRUE     # set FALSE only if you deliberately want v1.0.0
LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Resolve a table, preferring V1.1.0 over the root .csv.gz
pick <- function(stem) {
  v11  <- file.path(pic_dir, "V1.1.0", paste0(stem, ".csv"))
  upd  <- file.path(pic_dir, "update_to_v1.1.0", paste0(stem, ".csv"))
  root <- file.path(pic_dir, paste0(stem, ".csv.gz"))
  cands <- if (USE_V110) c(v11, upd, root) else c(root, v11, upd)
  for (p in cands) if (file.exists(p) && !dir.exists(p)) {
    cat(sprintf("  %-16s -> %s\n", stem, sub(paste0(pic_dir, "/"), "", p))); return(p)
  }
  stop("Could not find ", stem, " in V1.1.0/ or at root of: ", pic_dir)
}

cat("=== Resolving input files (preferring V1.1.0) ===\n")
f_icu  <- pick("ICUSTAYS")
f_pat  <- pick("PATIENTS")
f_adm  <- pick("ADMISSIONS")
f_lab  <- pick("LABEVENTS")
f_dlab <- pick("D_LABITEMS")

# ------------------------------------------------- 0. identify creatinine item
cat("\n=== Identifying the creatinine laboratory item ===\n")
dlab <- fread(f_dlab); setnames(dlab, tolower(names(dlab)))
lbl  <- intersect(c("label","labname","itemname","name"), names(dlab))[1]
cand <- dlab[grepl("creatinin|\u808c\u9150", get(lbl), ignore.case = TRUE)]
print(cand[, c("itemid", lbl), with = FALSE])
cand <- cand[!grepl("urin|\u5c3f|clearance|ratio", get(lbl), ignore.case = TRUE)]
CREAT_IDS <- unique(cand$itemid)
cat("Using itemid(s):", paste(CREAT_IDS, collapse = ", "), "\n")
if (!length(CREAT_IDS))
  stop("No serum creatinine itemid found. Inspect D_LABITEMS and set CREAT_IDS manually.")

# ---------------------------------------------------------------- 1. episodes
cat("\n=== Building episode table ===\n")
icu <- fread(f_icu); setnames(icu, tolower(names(icu)))
pat <- fread(f_pat); setnames(pat, tolower(names(pat)))
adm <- fread(f_adm); setnames(adm, tolower(names(adm)))
cat(sprintf("  ICUSTAYS rows read: %s   <-- confirm this against the release notes\n",
            format(nrow(icu), big.mark = ",")))

gv <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else rep(NA, nrow(dt))
icu[, `:=`(episode_id = as.character(icustay_id),
           patient_id = as.character(subject_id),
           hadm_id    = as.character(hadm_id),
           t0         = as.POSIXct(intime,  tz = "UTC"),
           icu_out    = as.POSIXct(outtime, tz = "UTC"))]
pat[, patient_id := as.character(subject_id)]
adm[, hadm_id    := as.character(hadm_id)]

ep <- merge(icu, pat, by = "patient_id", all.x = TRUE)
ep <- merge(ep,  adm, by = "hadm_id",    all.x = TRUE)

admt <- as.POSIXct(gv(ep, "admittime"), tz = "UTC")
dist <- as.POSIXct(gv(ep, "dischtime"), tz = "UTC")
dth  <- as.POSIXct(gv(ep, "deathtime"), tz = "UTC")
dob  <- as.POSIXct(gv(ep, "dob"),       tz = "UTC")
dod  <- suppressWarnings(as.POSIXct(gv(ep, "dod"), tz = "UTC"))

episodes <- data.table(
  source        = "PIC",
  episode_id    = ep$episode_id,
  patient_id    = ep$patient_id,
  hadm_id       = ep$hadm_id,
  t0            = ep$t0,
  # PIC stores date of birth; age in years at ICU admission
  age           = as.numeric(difftime(ep$t0, dob, units = "days")) / 365.25,
  sex           = as.character(gv(ep, "gender")),
  first_careunit= as.character(gv(ep, "first_careunit")),
  icu_los_days  = as.numeric(difftime(ep$icu_out, ep$t0, units = "days")),
  hosp_los_days = as.numeric(difftime(dist, admt, units = "days")),
  discharge_t_days = as.numeric(difftime(dist, ep$t0, units = "days")),
  died_inhosp   = as.integer(suppressWarnings(as.numeric(gv(ep, "hospital_expire_flag"))) == 1)
)
# in-hospital death time, falling back to date of death where absent.
# Both vectors are computed first so the assignment is fully vectorised
# (assigning a full-length RHS to a subset of rows is a data.table error).
dth_days <- as.numeric(difftime(dth, ep$t0, units = "days"))
dod_days <- as.numeric(difftime(dod, ep$t0, units = "days"))
episodes[, death_t_days := fifelse(is.na(dth_days), dod_days, dth_days)]
died_flag <- as.integer(!is.na(dth))
episodes[is.na(died_inhosp), died_inhosp := died_flag[is.na(died_inhosp)]]
episodes[age < 0 | age > 19, age := NA]

fwrite(episodes, paste0(out_dir, "PIC_episode_covariates.csv"))
cat(sprintf("  episodes written: %s (patients %s)\n",
            format(nrow(episodes), big.mark = ","),
            format(uniqueN(episodes$patient_id), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
cat("\n=== Extracting creatinine series ===\n")
con <- dbConnect(duckdb(), dbdir = ":memory:")
q <- sprintf("
  SELECT * FROM read_csv_auto('%s', header=true, strict_mode=false,
                              ignore_errors=true, quote='\"', escape='\"')
  WHERE ITEMID IN (%s)
", gsub("\\\\", "/", f_lab), paste(CREAT_IDS, collapse = ","))
cr <- as.data.table(dbGetQuery(con, q)); dbDisconnect(con, shutdown = TRUE)
setnames(cr, tolower(names(cr)))
cat(sprintf("  raw creatinine rows: %s\n", format(nrow(cr), big.mark = ",")))
if (!nrow(cr)) stop("No rows returned. Check CREAT_IDS against D_LABITEMS.")

v_col <- intersect(c("valuenum","value"), names(cr))[1]
t_col <- intersect(c("charttime","chart_time","storetime"), names(cr))[1]
s_col <- intersect(c("subject_id","patient_id"), names(cr))[1]
u_col <- intersect(c("valueuom","unit","units"), names(cr))[1]
cat(sprintf("  columns: subject=%s time=%s value=%s unit=%s\n", s_col, t_col, v_col, u_col))

cr <- cr[, .(patient_id = as.character(get(s_col)),
             charttime  = as.POSIXct(get(t_col), tz = "UTC"),
             value      = suppressWarnings(as.numeric(get(v_col))),
             unit_raw   = if (!is.na(u_col)) as.character(get(u_col)) else NA_character_)]
cr <- cr[!is.na(value) & !is.na(charttime) & value > 0]

cat("  units present in the source:\n"); print(table(cr$unit_raw, useNA = "ifany"))
med_raw <- median(cr$value)
cat(sprintf("  median raw value: %.2f\n", med_raw))
if (med_raw > 20) {
  cat("  -> umol/L detected; converting to mg/dL (divide by 88.42)\n")
  cr[, value := value / 88.42]
} else {
  cat("  -> mg/dL; no conversion applied\n")
}
# paediatric creatinine is LOW: neonatal values around 0.2-0.4 mg/dL are normal
cr <- cr[value > 0.02 & value < 30]

ser <- merge(cr, episodes[, .(episode_id, patient_id, t0)],
             by = "patient_id", allow.cartesian = TRUE)
ser[, t_days := as.numeric(difftime(charttime, t0, units = "days"))]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

series <- ser[, .(source = "PIC", episode_id, patient_id, charttime,
                  t_days = round(t_days, 5), value = round(value, 4), unit = "mg/dL")]
setorder(series, episode_id, t_days)
fwrite(series, paste0(out_dir, "PIC_creatinine_series.csv"))
cat(sprintf("  series written: %s rows\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
cat("\n=== Sanity checks (paste this into your reply) ===\n")
cat(sprintf("  median creatinine: %.2f mg/dL  (expect ~0.3-0.6 in children)\n", median(series$value)))
cat(sprintf("  age: median %.1f y (IQR %.1f-%.1f), max %.1f\n",
    median(episodes$age, na.rm = TRUE),
    quantile(episodes$age, .25, na.rm = TRUE),
    quantile(episodes$age, .75, na.rm = TRUE),
    max(episodes$age, na.rm = TRUE)))
cat(sprintf("  episodes with >=1 result: %s of %s (%.1f%%)\n",
    format(uniqueN(series$episode_id), big.mark = ","),
    format(nrow(episodes), big.mark = ","),
    100*uniqueN(series$episode_id)/nrow(episodes)))
n <- series[, .N, by = episode_id]
cat(sprintf("  results per episode: median %d (IQR %d-%d), max %d\n",
    as.integer(median(n$N)), as.integer(quantile(n$N,.25)),
    as.integer(quantile(n$N,.75)), max(n$N)))
cat(sprintf("  timing: pre-t0 %.1f%% | 7-365d before %.1f%% | 0-7d after %.1f%% | 7-90d after %.1f%%\n",
    100*mean(series$t_days < 0),
    100*mean(series$t_days <= -7 & series$t_days >= -365),
    100*mean(series$t_days > 0 & series$t_days <= 7),
    100*mean(series$t_days > 7 & series$t_days <= 90)))
cat(sprintf("  covariates: age %.1f%% | sex %.1f%% | hosp LOS %.1f%% | death time %.1f%%\n",
    100*mean(!is.na(episodes$age)), 100*mean(!is.na(episodes$sex)),
    100*mean(!is.na(episodes$hosp_los_days)), 100*mean(!is.na(episodes$death_t_days))))
cat(sprintf("  in-hospital mortality: %.1f%%\n", 100*mean(episodes$died_inhosp, na.rm = TRUE)))
cat(sprintf("  patients %s | episodes %s | episodes per patient %.2f\n",
    format(uniqueN(episodes$patient_id), big.mark = ","),
    format(nrow(episodes), big.mark = ","),
    nrow(episodes)/uniqueN(episodes$patient_id)))
cat("\n  CHECK: episode count should be close to the ICUSTAYS row count printed above.\n")
cat("         If it is much lower, a filter has been applied somewhere.\n")
cat("\nDone. Upload BOTH files. Do not open them in Excel.\n")
