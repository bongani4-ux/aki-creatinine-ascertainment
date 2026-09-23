

rm(list = ls())
library(data.table)
suppressWarnings(suppressMessages({ library(duckdb); library(DBI) }))

# ---------------------------------------------------------------- parameters
nwicu_dir <- "path/to/northwestern-icu-nwicu-database-0.1.0/northwestern-icu-nwicu-database-0.1.0"
out_dir   <- "path/to/export_series/"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# find a table whatever its extension / capitalisation / nesting
find_tbl <- function(dir, stem) {
  pat <- paste0("^", stem, "\\.csv(\\.gz)?$")
  hits <- list.files(dir, pattern = pat, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!dir.exists(hits)]
  if (!length(hits)) stop("Could not find '", stem, "' under: ", dir)
  cat(sprintf("  %-14s -> %s\n", stem, hits[1])); hits[1]
}

cat("=== Resolving input files ===\n")
f_icu  <- find_tbl(nwicu_dir, "icustays")
f_pat  <- find_tbl(nwicu_dir, "patients")
f_adm  <- find_tbl(nwicu_dir, "admissions")
f_lab  <- find_tbl(nwicu_dir, "labevents")
f_dlab <- find_tbl(nwicu_dir, "d_labitems")

# ------------------------------------------------- 0. identify creatinine item
cat("\n=== Identifying the creatinine lab item ===\n")
dlab <- fread(f_dlab); setnames(dlab, tolower(names(dlab)))
lbl  <- intersect(c("label","itemname","name"), names(dlab))[1]
cand <- dlab[grepl("creatinin", get(lbl), ignore.case = TRUE)]
print(cand)
# Serum creatinine only. Exclude urine creatinine and clearance measures.
cand <- cand[!grepl("urine|clearance|ratio|body fluid", get(lbl), ignore.case = TRUE)]
CREAT_IDS <- unique(cand$itemid)
cat("Using itemid(s):", paste(CREAT_IDS, collapse = ", "), "\n")
if (!length(CREAT_IDS)) stop("No serum creatinine itemid found. Inspect d_labitems and set CREAT_IDS manually.")

# ---------------------------------------------------------------- 1. episodes
cat("\n=== Building episode table ===\n")
icu <- fread(f_icu); setnames(icu, tolower(names(icu)))
pat <- fread(f_pat); setnames(pat, tolower(names(pat)))
adm <- fread(f_adm); setnames(adm, tolower(names(adm)))

idc <- function(dt, opts) intersect(opts, names(dt))[1]
c_stay <- idc(icu, c("stay_id","icustay_id"))
c_subj <- idc(icu, c("subject_id","patient_id"))
c_hadm <- idc(icu, c("hadm_id","visit_id"))

icu[, `:=`(episode_id = as.character(get(c_stay)),
           patient_id = as.character(get(c_subj)),
           hadm_id    = as.character(get(c_hadm)),
           t0         = as.POSIXct(intime,  tz = "UTC"),
           icu_out    = as.POSIXct(outtime, tz = "UTC"))]

pat[, patient_id := as.character(get(idc(pat, c("subject_id","patient_id"))))]
adm[, hadm_id    := as.character(get(idc(adm, c("hadm_id","visit_id"))))]

ep <- merge(icu, pat, by = "patient_id", all.x = TRUE)
ep <- merge(ep,  adm, by = "hadm_id",    all.x = TRUE)

gv <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else rep(NA, nrow(dt))
adt <- as.POSIXct(gv(ep, "admittime"), tz = "UTC")
ddt <- as.POSIXct(gv(ep, "dischtime"), tz = "UTC")
dth <- as.POSIXct(gv(ep, "deathtime"), tz = "UTC")
dod <- suppressWarnings(as.POSIXct(gv(ep, "dod"), tz = "UTC"))

episodes <- data.table(
  source        = "NWICU",
  episode_id    = ep$episode_id,
  patient_id    = ep$patient_id,
  hadm_id       = ep$hadm_id,
  t0            = ep$t0,
  age           = suppressWarnings(as.numeric(gv(ep, "anchor_age"))),
  sex           = as.character(gv(ep, "gender")),
  race          = as.character(gv(ep, "race")),
  admission_type= as.character(gv(ep, "admission_type")),
  first_careunit= as.character(gv(ep, "first_careunit")),
  icu_los_days  = as.numeric(difftime(ep$icu_out, ep$t0, units = "days")),
  hosp_los_days = as.numeric(difftime(ddt, adt, units = "days")),
  discharge_t_days = as.numeric(difftime(ddt, ep$t0, units = "days")),
  died_inhosp   = as.integer(suppressWarnings(as.numeric(gv(ep, "hospital_expire_flag"))) == 1)
)
episodes[, death_t_days := as.numeric(difftime(dth, ep$t0, units = "days"))]
episodes[is.na(death_t_days) & !is.na(dod),
         death_t_days := as.numeric(difftime(dod, ep$t0, units = "days"))]
episodes[is.na(died_inhosp), died_inhosp := as.integer(!is.na(dth))]
if (all(is.na(episodes$age)) && "age_years" %in% names(ep))
  episodes[, age := suppressWarnings(as.numeric(ep$age_years))]

fwrite(episodes, paste0(out_dir, "NWICU_episode_covariates.csv"))
cat(sprintf("  episodes written: %s\n", format(nrow(episodes), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
cat("\n=== Extracting creatinine series ===\n")
con <- dbConnect(duckdb(), dbdir = ":memory:")
q <- sprintf("
  SELECT * FROM read_csv_auto('%s', header=true, strict_mode=false,
                              ignore_errors=true, quote='\"', escape='\"')
  WHERE itemid IN (%s)
", gsub("\\\\", "/", f_lab), paste(CREAT_IDS, collapse = ","))
cr <- as.data.table(dbGetQuery(con, q)); dbDisconnect(con, shutdown = TRUE)
setnames(cr, tolower(names(cr)))
cat(sprintf("  raw creatinine rows: %s\n", format(nrow(cr), big.mark = ",")))

v_col <- intersect(c("valuenum","value"), names(cr))[1]
t_col <- intersect(c("charttime","storetime","chart_time"), names(cr))[1]
s_col <- intersect(c("subject_id","patient_id"), names(cr))[1]
u_col <- intersect(c("valueuom","unit","units"), names(cr))[1]
cat(sprintf("  columns: subject=%s time=%s value=%s unit=%s\n", s_col, t_col, v_col, u_col))

cr <- cr[, .(patient_id = as.character(get(s_col)),
             charttime  = as.POSIXct(get(t_col), tz = "UTC"),
             value      = suppressWarnings(as.numeric(get(v_col))),
             unit_raw   = if (!is.na(u_col)) as.character(get(u_col)) else NA_character_)]
cr <- cr[!is.na(value) & !is.na(charttime) & value > 0]

# ---- UNIT HANDLING: declare, do not assume ----
cat("  units present in the source:\n"); print(table(cr$unit_raw, useNA = "ifany"))
med_raw <- median(cr$value)
cat(sprintf("  median raw value: %.2f\n", med_raw))
if (med_raw > 20) {
  cat("  -> looks like umol/L; converting to mg/dL (divide by 88.42)\n")
  cr[, value := value / 88.42]
} else {
  cat("  -> looks like mg/dL; no conversion applied\n")
}
cr <- cr[value > 0.05 & value < 50]

ser <- merge(cr, episodes[, .(episode_id, patient_id, t0)],
             by = "patient_id", allow.cartesian = TRUE)
ser[, t_days := as.numeric(difftime(charttime, t0, units = "days"))]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

series <- ser[, .(source = "NWICU", episode_id, patient_id, charttime,
                  t_days = round(t_days, 5), value = round(value, 4), unit = "mg/dL")]
setorder(series, episode_id, t_days)
fwrite(series, paste0(out_dir, "NWICU_creatinine_series.csv"))
cat(sprintf("  series written: %s rows\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
source_checks <- function(series, episodes) {
  cat("\n=== Sanity checks (paste this into your reply) ===\n")
  cat(sprintf("  median creatinine: %.2f mg/dL  (expect ~0.7-1.2)\n", median(series$value)))
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
  cat(sprintf("  patients %s | episodes %s | episodes per patient %.2f\n",
      format(uniqueN(episodes$patient_id), big.mark = ","),
      format(nrow(episodes), big.mark = ","),
      nrow(episodes)/uniqueN(episodes$patient_id)))
}
source_checks(series, episodes)
cat("\nDone. Upload BOTH files. Do not open them in Excel.\n")
