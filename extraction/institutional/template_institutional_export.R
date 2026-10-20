# =============================================================================
# Institutional surgical cohorts (ZD, NJF, AHMU) - specification template
#
# This is a TEMPLATE for the hospital data team. It assumes the raw operation
# and laboratory tables rather than the previously supplied summary extract.
#
# Produces, per cohort:
#     <COHORT>_creatinine_series.csv
#     <COHORT>_episode_covariates.csv
#
# IMPORTANT — two defects in the previous extract that this export must fix:
#   1. ZD PATIENT_ID was truncated to the literal string 'ZY', and VISIT_ID and
#      OPER_ID were null in all 312,579 rows. Identifiers must be exported
#      intact (as text, not numeric) so episodes can be linked and deduplicated.
#   2. END_PLACE was empty in every row of both cohorts, so mortality could not
#      be determined. A death date or in-hospital death flag is required.
# =============================================================================

rm(list = ls())
library(data.table)

# ---------------------------------------------------------------- parameters
COHORT  <- "ZD"          # set to "ZD" or "NJF" and run twice
base_dir <- "path/to/de-identified operations file/"
out_dir  <- "path/to/export_series/"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90

# Source tables. 01 = demographics, 03 = operation detail, 13 = ward/ICU times.
# The laboratory table is NOT among the files supplied so far and must be added.
f_demo <- file.path(base_dir, COHORT, "01.csv")
f_oper <- file.path(base_dir, COHORT, "03.csv")
f_ward <- file.path(base_dir, COHORT, "13.csv")
f_labs <- file.path(base_dir, COHORT, "LABS.csv")   # <-- TO BE SUPPLIED

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
KEY <- c("PATIENT_ID", "VISIT_ID", "OPER_ID")

# helper: strip the stray single quotes seen in the supplied CSVs
unq <- function(x) gsub("^'|'$", "", trimws(as.character(x)))

# ---------------------------------------------------------------- 1. episodes
cat("=== Building episode table for", COHORT, "===\n")
demo <- fread(f_demo); oper <- fread(f_oper); ward <- fread(f_ward)
for (dt in list(demo, oper, ward)) for (k in KEY) if (k %in% names(dt)) dt[, (k) := unq(get(k))]

demo <- unique(demo, by = KEY); oper <- unique(oper, by = KEY); ward <- unique(ward, by = KEY)

ep <- merge(oper, demo, by = KEY, all.x = TRUE)
ep <- merge(ep,   ward, by = KEY, all.x = TRUE)

# ---- time zero: anaesthesia/incision start (C1) ----
ep[, t0 := as.POSIXct(START_DATE_TIME, tz = "UTC")]
ep[is.na(t0), t0 := as.POSIXct(ANES_START_TIME, tz = "UTC")]
ep[is.na(t0), t0 := as.POSIXct(IN_DATE_TIME,   tz = "UTC")]

# ---- sex: the raw field mixes Chinese, codes and stray values ----
clean_sex <- function(x) {
  s <- trimws(as.character(x)); out <- rep(NA_character_, length(s))
  out[grepl("\u7537", s)] <- "Male"      # 男
  out[grepl("\u5973", s)] <- "Female"    # 女
  out[s %in% c("1", "M", "Male")]   <- "Male"
  out[s %in% c("2", "F", "Female")] <- "Female"
  out
}
# ---- ASA: mixed Roman, Chinese numerals and a trailing 级 character ----
clean_asa <- function(x) {
  s <- gsub("\u7ea7", "", trimws(as.character(x)))     # strip 级
  out <- rep(NA_integer_, length(s))
  out[s %in% c("I","\u2160","\u2170","1","\uff11")] <- 1L
  out[s %in% c("II","\u2161","\u2171","2","\uff12")] <- 2L
  out[s %in% c("III","\u2162","\u2172","3","\uff13")] <- 3L
  out[s %in% c("IV","\u2163","\u2173","4","\uff14")] <- 4L
  out[s %in% c("V","\u2164","\u2174","5","\uff15")] <- 5L
  out
}

hosp_in  <- as.POSIXct(ep$IN_HOSPITAL_TIME,  tz = "UTC")
hosp_out <- as.POSIXct(ep$OUT_HOSPITAL_TIME, tz = "UTC")
icu_in   <- as.POSIXct(ep$IN_ICU_TIME,       tz = "UTC")
icu_out  <- as.POSIXct(ep$OUT_ICU_TIME,      tz = "UTC")

# ---- mortality: REQUIRED, currently missing ----
# END_PLACE was empty in the supplied files. Supply ONE of:
#   DEATH_DATE_TIME   a date-time of death, or
#   INHOSP_DEATH_FLAG 0/1, or
#   END_PLACE         populated discharge-outcome code with a death category
death_time <- if ("DEATH_DATE_TIME" %in% names(ep))
                as.POSIXct(ep$DEATH_DATE_TIME, tz = "UTC") else as.POSIXct(NA)
died_flag  <- if ("INHOSP_DEATH_FLAG" %in% names(ep)) as.integer(ep$INHOSP_DEATH_FLAG) else NA_integer_

episodes <- data.table(
  source        = COHORT,
  episode_id    = paste(ep$PATIENT_ID, ep$VISIT_ID, ep$OPER_ID, sep = "|"),
  patient_id    = ep$PATIENT_ID,
  visit_id      = ep$VISIT_ID,
  oper_id       = ep$OPER_ID,
  t0            = ep$t0,
  age           = suppressWarnings(as.numeric(ep$AGE)),
  sex           = clean_sex(ep$SEX),
  bmi           = suppressWarnings(as.numeric(ep$BMI)),
  asa           = clean_asa(ep$ASA_GRADE),
  emergency     = suppressWarnings(as.numeric(ep$EMERGENCY_INDICATOR)),
  operation_scale = as.character(ep$OPERATION_SCALE),
  department    = as.character(ep$OPERATING_DEPT),
  hosp_los_days = as.numeric(difftime(hosp_out, hosp_in, units = "days")),
  icu_los_days  = as.numeric(difftime(icu_out,  icu_in,  units = "days")),
  icu_admitted  = as.integer(!is.na(icu_in)),
  discharge_t_days = as.numeric(difftime(hosp_out, ep$t0, units = "days")),
  death_t_days     = as.numeric(difftime(death_time, ep$t0, units = "days")),
  died_inhosp      = died_flag
)
episodes[age < 0 | age > 110, age := NA]
episodes[hosp_los_days < 0 | hosp_los_days > 365, hosp_los_days := NA]

fwrite(episodes, paste0(out_dir, COHORT, "_episode_covariates.csv"))
cat(sprintf("  episodes: %s rows\n", format(nrow(episodes), big.mark = ",")))
cat(sprintf("  identifiers intact: %.1f%%   age: %.1f%%   sex: %.1f%%   ASA: %.1f%%   death: %.1f%%\n",
  100*mean(nchar(episodes$patient_id) > 3),
  100*mean(!is.na(episodes$age)), 100*mean(!is.na(episodes$sex)),
  100*mean(!is.na(episodes$asa)), 100*mean(!is.na(episodes$death_t_days))))
if (mean(nchar(episodes$patient_id) > 3) < 0.9)
  warning("PATIENT_ID still looks truncated. Export identifiers as TEXT.")

# ---------------------------------------------------------------- 2. series
cat("\n=== Extracting creatinine series ===\n")
if (!file.exists(f_labs)) {
  cat("  LABS.csv not found.\n")
  cat("  Ask the data team for a laboratory table containing, for every\n")
  cat("  creatinine result on these patients between 365 days before and\n")
  cat("  90 days after each operation:\n")
  cat("      PATIENT_ID, VISIT_ID (if available), TEST_DATE_TIME,\n")
  cat("      ITEM_NAME (or code), RESULT_VALUE, RESULT_UNIT\n")
  cat("  Outpatient and readmission results must be included, not just\n")
  cat("  those from the index admission.\n")
  quit(save = "no")
}

labs <- fread(f_labs)
for (k in intersect(KEY, names(labs))) labs[, (k) := unq(get(k))]
setnames(labs, tolower(names(labs)))

item_col  <- intersect(c("item_name","test_name","itemname","item"), names(labs))[1]
time_col  <- intersect(c("test_date_time","result_time","charttime","test_time"), names(labs))[1]
value_col <- intersect(c("result_value","value","result"), names(labs))[1]
unit_col  <- intersect(c("result_unit","unit","units"), names(labs))[1]
cat(sprintf("  columns: item=%s time=%s value=%s unit=%s\n",
            item_col, time_col, value_col, unit_col))

# Chinese creatinine labels: 肌酐 / 血肌酐 / 血清肌酐
cr <- labs[grepl("\u808c\u9150|creat", get(item_col), ignore.case = TRUE),
           .(patient_id = get("patient_id"),
             chart_time = as.POSIXct(get(time_col), tz = "UTC"),
             value      = suppressWarnings(as.numeric(get(value_col))),
             unit       = if (!is.na(unit_col)) as.character(get(unit_col)) else NA_character_)]
cr <- cr[!is.na(value) & value > 0]

# harmonise units: institutional labs report umol/L (C13)
cr[, unit_std := fifelse(grepl("mg", unit, ignore.case = TRUE), "mg/dL", "umol/L")]
cr[, value_umol := fifelse(unit_std == "mg/dL", value * 88.42, value)]
cr <- cr[value_umol > 5 & value_umol < 4000]

ser <- merge(cr, episodes[, .(episode_id, patient_id, t0)],
             by = "patient_id", allow.cartesian = TRUE)
ser[, t_days := as.numeric(difftime(chart_time, t0, units = "days"))]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

# The analysis pipeline works in mg/dL (docs/data_specification.md), so the
# harmonised umol/L value is converted back on output.
series <- ser[, .(source = COHORT, episode_id, patient_id, chart_time,
                  t_days = round(t_days, 5),
                  value = round(value_umol / 88.42, 4), unit = "mg/dL")]
setorder(series, episode_id, t_days)
fwrite(series, paste0(out_dir, COHORT, "_creatinine_series.csv"))
cat(sprintf("  series: %s rows\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
cat("\n=== Sanity checks ===\n")
cat(sprintf("  median creatinine: %.2f mg/dL (expect ~0.7-1.0)\n", median(series$value)))
n_per <- series[, .N, by = episode_id]
cat(sprintf("  episodes with >=1 result: %s (%.1f%%)\n",
    format(uniqueN(series$episode_id), big.mark = ","),
    100*uniqueN(series$episode_id)/nrow(episodes)))
cat(sprintf("  results per episode: median %d (IQR %d-%d)\n",
    as.integer(median(n_per$N)), as.integer(quantile(n_per$N,.25)), as.integer(quantile(n_per$N,.75))))
cat(sprintf("  pre-t0: %.1f%%   7-30d before t0: %.1f%%   day 7-90 after: %.1f%%\n",
    100*mean(series$t_days < 0),
    100*mean(series$t_days <= -7 & series$t_days >= -30),
    100*mean(series$t_days > 7 & series$t_days <= 90)))

cat("\nDone. Upload BOTH files for this cohort, then set COHORT to the other and re-run.\n")
