# =============================================================================
# EXPORT — Zigong  (Zigong Fourth People's Hospital infection cohort, PhysioNet)
# Full creatinine measurement series + episode-level covariates
#
# FILES NEEDED (all inside DataTables/):
#     dtBaseline.csv      one row per admission: demographics, admit/discharge times
#     dtLab.csv           all laboratory results  (150 MB uncompressed)
#     dtOutCome.csv       discharge outcome, including death
#     datDictionary.csv   variable dictionary (used to confirm the lab item name)
#
# Produces:
#     Zigong_creatinine_series.csv
#     Zigong_episode_covariates.csv
#
# t0 = hospital admission.
#
# NOTE: this cohort is restricted to patients admitted WITH INFECTION, so it is
# not a general ICU population. Keep that in mind when interpreting its
# ascertainment relative to the other sources.
# =============================================================================

rm(list = ls())
library(data.table)

# ---------------------------------------------------------------- parameters
zigong_dir <- "path/to/Critical care database comprising patients with infection at Zigong Fourth People's Hospital"
out_dir    <- "path/to/export_series/"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

find_one <- function(dir, stem) {
  h <- list.files(dir, pattern = paste0("^", stem, "\\.csv$"),
                  recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  h <- h[!dir.exists(h)]
  if (!length(h)) stop("Could not find '", stem, ".csv' under: ", dir)
  cat(sprintf("  %-16s -> %s\n", stem, h[1])); h[1]
}

cat("=== Resolving input files ===\n")
f_base <- find_one(zigong_dir, "dtBaseline")
f_lab  <- find_one(zigong_dir, "dtLab")
f_out  <- find_one(zigong_dir, "dtOutCome")

# ------------------------------------------------- 0. identify creatinine item
cat("\n=== Identifying the creatinine laboratory item ===\n")
lab_head <- fread(f_lab, nrows = 200000)
cat("dtLab columns:\n"); print(names(lab_head))
setnames(lab_head, tolower(names(lab_head)))
item_col_h <- intersect(c("item_name","itemname","labname","item","test_name"),
                        names(lab_head))[1]
if (is.na(item_col_h)) stop("Could not find an item-name column in dtLab. Columns are:\n  ",
                            paste(names(lab_head), collapse = ", "))
cat(sprintf("\nDistinct labels matching creatinine (column '%s'):\n", item_col_h))
# Chinese: 肌酐 / 血肌酐 / 血清肌酐 ; English: creatinine
hits <- unique(grep("\u808c\u9150|creatinin", lab_head[[item_col_h]],
                    ignore.case = TRUE, value = TRUE))
print(hits)
if (!length(hits)) {
  cat("\n!! Nothing matched. First 60 distinct labels in dtLab:\n")
  print(head(unique(lab_head[[item_col_h]]), 60))
  stop("Set CREAT_PATTERN below to the correct label and re-run.")
}
# exclude urine creatinine and clearance if present
CREAT_PATTERN <- "\u808c\u9150|creatinin"
EXCLUDE_PATTERN <- "\u5c3f|urin|clearance|ratio|\u6e05\u9664"

# ---------------------------------------------------------------- 1. episodes
cat("\n=== Building episode table ===\n")
base <- fread(f_base); setnames(base, tolower(names(base)))
outc <- fread(f_out);  setnames(outc, tolower(names(outc)))
cat("dtBaseline columns:\n"); print(names(base))
cat("dtOutCome columns:\n");  print(names(outc))

gv <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else rep(NA, nrow(dt))
key <- intersect(c("inp_no","patient_id","id"), names(base))[1]
cat(sprintf("  joining on: %s\n", key))

if (key %in% names(outc)) {
  base[[key]] <- as.character(base[[key]])
  outc[[key]] <- as.character(outc[[key]])
  base <- merge(base, outc, by = key, all.x = TRUE, suffixes = c("", "_out"))
}

adm <- as.POSIXct(gv(base, "admission_time"), tz = "UTC")
if (all(is.na(adm))) adm <- as.POSIXct(gv(base, "admittime"), tz = "UTC")
if (all(is.na(adm))) adm <- as.POSIXct(gv(base, "in_hospital_time"), tz = "UTC")
dis <- as.POSIXct(gv(base, "discharge_time"), tz = "UTC")
if (all(is.na(dis))) dis <- as.POSIXct(gv(base, "discharge_date_time"), tz = "UTC")

# death: dtOutCome carries the discharge outcome
out_raw <- as.character(gv(base, "outcome"))
if (all(is.na(out_raw))) out_raw <- as.character(gv(base, "discharge_status"))
died <- as.integer(grepl("death|died|dead|\u6b7b\u4ea1", out_raw, ignore.case = TRUE))

episodes <- data.table(
  source        = "Zigong",
  episode_id    = as.character(base[[key]]),
  patient_id    = as.character(gv(base, "patient_id")),
  t0            = 0,   # LabTime offsets are measured from admission
  admission_time = adm,
  age           = suppressWarnings(as.numeric(gv(base, "age"))),
  sex           = as.character(gv(base, "sex")),
  admit_dept    = as.character(gv(base, "admitdept")),
  infection_site= as.character(gv(base, "infectionsite")),
  discharge_outcome = out_raw,
  died_inhosp   = died
)
if (all(is.na(episodes$patient_id))) episodes[, patient_id := episode_id]
episodes[, `:=`(
  hosp_los_days    = as.numeric(difftime(dis, adm, units = "days")),
  discharge_t_days = as.numeric(difftime(dis, adm, units = "days")),
  icu_los_days     = NA_real_,
  # Zigong records the outcome but not a death timestamp; approximate the time
  # of death by the discharge time for those who died, and leave it NA otherwise.
  death_t_days     = ifelse(died == 1, as.numeric(difftime(dis, adm, units = "days")), NA_real_)
)]
episodes[age < 0 | age > 110, age := NA]
episodes[hosp_los_days < 0 | hosp_los_days > 365, hosp_los_days := NA]

fwrite(episodes, paste0(out_dir, "Zigong_episode_covariates.csv"))
cat(sprintf("  episodes written: %s\n", format(nrow(episodes), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
# NOTE ON TIME: dtLab$LabTime is NOT a timestamp. It is a numeric offset in
# DAYS from hospital admission (e.g. 0.600000 = 14.4 h after admission).
# So t_days = LabTime directly; no datetime arithmetic is required.
cat("\n=== Extracting creatinine series ===\n")
lab <- fread(f_lab); setnames(lab, tolower(names(lab)))
if ("v1" %in% names(lab)) lab[, v1 := NULL]          # drop the row-index column
cat("  dtLab columns: "); cat(paste(names(lab), collapse = ", "), "\n")

item_col <- intersect(c("item","item_name","itemname","labname"), names(lab))[1]
time_col <- intersect(c("labtime","time","record_time","charttime"), names(lab))[1]
val_col  <- intersect(c("labvalue","value","result"), names(lab))[1]
unit_col <- intersect(c("unit","units"), names(lab))[1]
key_col  <- intersect(c("inp_no","patient_id","id"), names(lab))[1]
cat(sprintf("  columns: key=%s item=%s time=%s value=%s unit=%s\n",
            key_col, item_col, time_col, val_col, unit_col))
if (any(is.na(c(item_col, time_col, val_col, key_col))))
  stop("Could not identify a required dtLab column. Names are:\n  ",
       paste(names(lab), collapse = ", "))

cr <- lab[grepl(CREAT_PATTERN, get(item_col), ignore.case = TRUE) &
          !grepl(EXCLUDE_PATTERN, get(item_col), ignore.case = TRUE),
          .(episode_id = as.character(get(key_col)),
            t_days     = suppressWarnings(as.numeric(get(time_col))),
            value      = suppressWarnings(as.numeric(get(val_col))),
            unit_raw   = as.character(get(unit_col)))]
cr <- cr[!is.na(value) & !is.na(t_days) & value > 0]
cat(sprintf("  raw creatinine rows: %s\n", format(nrow(cr), big.mark = ",")))
if (!nrow(cr)) stop("No creatinine rows found. Check CREAT_PATTERN against the printed labels.")

cat("  item labels captured:\n")
print(unique(lab[grepl(CREAT_PATTERN, get(item_col), ignore.case = TRUE) &
                 !grepl(EXCLUDE_PATTERN, get(item_col), ignore.case = TRUE)][[item_col]]))
cat("  units present:\n"); print(table(cr$unit_raw, useNA = "ifany"))
cat(sprintf("  LabTime range (days from admission): %.2f to %.2f\n",
            min(cr$t_days), max(cr$t_days)))

# ---- UNIT HANDLING: Chinese labs normally report umol/L ----
med_raw <- median(cr$value)
cat(sprintf("  median raw value: %.1f\n", med_raw))
if (med_raw > 20) {
  cat("  -> umol/L detected; converting to mg/dL (divide by 88.42)\n")
  cr[, value := value / 88.42]
} else {
  cat("  -> mg/dL; no conversion applied\n")
}
cr <- cr[value > 0.05 & value < 50]

cr <- cr[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]
cr <- merge(cr, episodes[, .(episode_id, patient_id)], by = "episode_id", all.x = TRUE)

series <- cr[, .(source = "Zigong", episode_id, patient_id,
                 charttime = NA_character_,     # Zigong supplies offsets only
                 t_days = round(t_days, 5),
                 value  = round(value, 4), unit = "mg/dL")]
setorder(series, episode_id, t_days)
fwrite(series, paste0(out_dir, "Zigong_creatinine_series.csv"))
cat(sprintf("  series written: %s rows\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
cat("\n=== Sanity checks (paste this into your reply) ===\n")
cat(sprintf("  median creatinine: %.2f mg/dL  (expect ~0.7-1.4 in an infection cohort)\n",
            median(series$value)))
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
cat(sprintf("  covariates: age %.1f%% | sex %.1f%% | hosp LOS %.1f%% | died flag %.1f%%\n",
    100*mean(!is.na(episodes$age)), 100*mean(!is.na(episodes$sex)),
    100*mean(!is.na(episodes$hosp_los_days)), 100*mean(!is.na(episodes$died_inhosp))))
cat(sprintf("  in-hospital mortality: %.1f%%\n", 100*mean(episodes$died_inhosp, na.rm = TRUE)))
cat(sprintf("  patients %s | episodes %s | episodes per patient %.2f\n",
    format(uniqueN(episodes$patient_id), big.mark = ","),
    format(nrow(episodes), big.mark = ","),
    nrow(episodes)/uniqueN(episodes$patient_id)))
cat("\n  NOTE: death time is approximated by discharge time for those who died,\n")
cat("        because Zigong records the outcome but not a death timestamp.\n")
cat("\nDone. Upload BOTH files. Do not open them in Excel.\n")
