# =============================================================================
# EXPORT — HiRID  (Bern high time-resolution ICU dataset, PhysioNet)
# Full creatinine measurement series + episode-level covariates
#
# FILES NEEDED:
#     general_table.csv                     (one row per admission)
#     hirid_variable_reference.csv          (maps variableid -> concept name)
#     observation_tables/  (parquet OR csv)  (the actual measurements)
#
# Produces:
#     HiRID_creatinine_series.csv
#     HiRID_episode_covariates.csv
#
# t0 = ICU admission (admissiontime in general_table)
#
# NOTE: HiRID contains ICU data only. It has no pre-admission laboratory
# history and ends at ICU discharge. Expect almost all results to fall between
# t0 and ICU discharge, with nothing before t0 and little after day 7.
# That is a property of the source, not an extraction error - export it anyway,
# because the emptiness of those windows is exactly what we are quantifying.
# =============================================================================

rm(list = ls())
library(data.table)
suppressWarnings(suppressMessages({ library(duckdb); library(DBI) }))

# ---------------------------------------------------------------- parameters
hirid_dir <- "path/to/hirid-a-high-time-resolution-icu-dataset-1.1.1/hirid-a-high-time-resolution-icu-dataset-1.1.1"
out_dir   <- "path/to/export_series/"

# HiRID ships the measurement tables as tar.gz archives. They must be unpacked
# before DuckDB can read them. This is where they will be extracted to.
# Needs roughly 10 GB free. Extraction runs once; later runs reuse it.
work_dir  <- "path/to/hirid_extracted/"

LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90

if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)

find_one <- function(dir, pat) {
  h <- list.files(dir, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  h <- h[!dir.exists(h)]
  if (!length(h)) stop("Could not find a file matching '", pat, "' under: ", dir)
  cat(sprintf("  %-30s -> %s\n", pat, h[1])); h[1]
}

cat("=== Resolving reference files ===\n")
# reference_data is already unpacked in the download; if not, unpack it first
if (!length(list.files(hirid_dir, pattern = "^general_table\\.csv$",
                       recursive = TRUE, ignore.case = TRUE))) {
  ref_tar <- find_one(hirid_dir, "^reference_data\\.tar\\.gz$")
  cat("  unpacking reference_data.tar.gz ...\n")
  untar(ref_tar, exdir = work_dir)
}
search_dirs <- c(hirid_dir, work_dir)
f_gen <- NULL; f_ref <- NULL
for (d in search_dirs) {
  if (is.null(f_gen)) try(f_gen <- find_one(d, "^general_table\\.csv(\\.gz)?$"), silent = TRUE)
  if (is.null(f_ref)) try(f_ref <- find_one(d, "hirid_variable_reference\\.csv$"), silent = TRUE)
}
if (is.null(f_gen) || is.null(f_ref))
  stop("general_table.csv or hirid_variable_reference.csv not found in either\n  ",
       hirid_dir, "\n  ", work_dir)

# ---------------------------------------------------------------- unpack observations
cat("\n=== Locating observation tables ===\n")
obs_root <- file.path(work_dir, "observation_tables")
pq <- list.files(obs_root, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
cs <- list.files(obs_root, pattern = "\\.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)

if (!length(pq) && !length(cs)) {
  # prefer parquet: smaller to query and typed correctly
  tar_pq <- list.files(hirid_dir, pattern = "observation_tables_parquet\\.tar\\.gz$",
                       recursive = TRUE, full.names = TRUE)
  tar_cs <- list.files(hirid_dir, pattern = "observation_tables_csv\\.tar\\.gz$",
                       recursive = TRUE, full.names = TRUE)
  tar_use <- if (length(tar_pq)) tar_pq[1] else if (length(tar_cs)) tar_cs[1] else NULL
  if (is.null(tar_use)) stop("No observation_tables archive found under: ", hirid_dir)
  cat(sprintf("  unpacking %s\n", basename(tar_use)))
  cat("  this is ~6.6 GB and takes several minutes; it only happens once\n")
  untar(tar_use, exdir = work_dir)
  pq <- list.files(work_dir, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  cs <- list.files(work_dir, pattern = "\\.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  cs <- cs[grepl("observation", cs, ignore.case = TRUE)]
}
pq <- pq[grepl("observation", pq, ignore.case = TRUE)]
cat(sprintf("  found %d parquet and %d csv part files\n", length(pq), length(cs)))
if (!length(pq) && !length(cs)) stop("Extraction produced no observation part files.")
obs_glob <- if (length(pq)) dirname(pq[1]) else dirname(cs[1])
use_parquet <- length(pq) > 0
cat(sprintf("  querying: %s (%s)\n", obs_glob, if (use_parquet) "parquet" else "csv"))

# ------------------------------------------------- 0. identify creatinine variable
cat("\n=== Identifying the creatinine variable ===\n")
ref <- fread(f_ref); setnames(ref, tolower(names(ref)))
nm  <- intersect(c("variablename","variable name","uid","name","concept"), names(ref))[1]
idc <- intersect(c("variableid","variable id","id"), names(ref))[1]
cand <- ref[grepl("creatinin|kreatinin", get(nm), ignore.case = TRUE)]
print(unique(cand[, c(idc, nm), with = FALSE]))
cand <- cand[!grepl("urin|clearance|ratio", get(nm), ignore.case = TRUE)]
CREAT_IDS <- unique(suppressWarnings(as.numeric(cand[[idc]])))
CREAT_IDS <- CREAT_IDS[!is.na(CREAT_IDS)]
cat("Using variableid(s):", paste(CREAT_IDS, collapse = ", "), "\n")
if (!length(CREAT_IDS))
  stop("No serum creatinine variableid found. Inspect the reference file and set CREAT_IDS manually (serum creatinine is commonly 20000600).")

# ---------------------------------------------------------------- 1. episodes
cat("\n=== Building episode table ===\n")
gen <- fread(f_gen); setnames(gen, tolower(names(gen)))
gv <- function(dt, nm) if (nm %in% names(dt)) dt[[nm]] else rep(NA, nrow(dt))

t0 <- as.POSIXct(gv(gen, "admissiontime"), tz = "UTC")
dstat <- as.character(gv(gen, "discharge_status"))

episodes <- data.table(
  source        = "HiRID",
  episode_id    = as.character(gen$patientid),
  patient_id    = as.character(gen$patientid),   # HiRID: one admission per patientid
  t0            = t0,
  age           = suppressWarnings(as.numeric(gv(gen, "age"))),
  sex           = as.character(gv(gen, "sex")),
  discharge_status = dstat,
  died_inhosp   = as.integer(tolower(dstat) %in% c("dead","deceased","died"))
)
# HiRID gives no discharge timestamp and no death timestamp.
# Leave the columns present but empty so the pooled analysis knows they are
# absent rather than zero (criterion C11 cannot be applied to this source).
episodes[, `:=`(hosp_los_days = NA_real_, icu_los_days = NA_real_,
                discharge_t_days = NA_real_, death_t_days = NA_real_)]

fwrite(episodes, paste0(out_dir, "HiRID_episode_covariates.csv"))
cat(sprintf("  episodes written: %s\n", format(nrow(episodes), big.mark = ",")))

# ---------------------------------------------------------------- 2. series
cat("\n=== Extracting creatinine series ===\n")
con <- dbConnect(duckdb(), dbdir = ":memory:")
ids <- paste(format(CREAT_IDS, scientific = FALSE), collapse = ",")
src <- if (use_parquet) {
  sprintf("read_parquet('%s/**/*.parquet', union_by_name=true)", gsub("\\\\", "/", obs_glob))
} else {
  sprintf("read_csv_auto('%s/**/*.csv*', union_by_name=true, header=true, ignore_errors=true)",
          gsub("\\\\", "/", obs_glob))
}
q <- sprintf("SELECT patientid, datetime, value, variableid FROM %s WHERE variableid IN (%s)", src, ids)
cr <- as.data.table(dbGetQuery(con, q)); dbDisconnect(con, shutdown = TRUE)
setnames(cr, tolower(names(cr)))
cat(sprintf("  raw creatinine rows: %s\n", format(nrow(cr), big.mark = ",")))
if (!nrow(cr)) stop("No rows returned. Check CREAT_IDS against the reference file.")

cr <- cr[, .(patient_id = as.character(patientid),
             charttime  = as.POSIXct(datetime, tz = "UTC"),
             value      = suppressWarnings(as.numeric(value)))]
cr <- cr[!is.na(value) & !is.na(charttime) & value > 0]
# The published export contained each draw 2-5 times (same patient, second and
# value: 114,181 of 190,520 rows), from overlapping part files under the glob
# and/or more than one matching variableid. Keep one row per draw.
n_before <- nrow(cr)
cr <- unique(cr, by = c("patient_id", "charttime", "value"))
cat(sprintf("  exact duplicate rows removed: %s\n", format(n_before - nrow(cr), big.mark = ",")))

# ---- UNIT HANDLING: HiRID reports creatinine in umol/L ----
med_raw <- median(cr$value)
cat(sprintf("  median raw value: %.1f\n", med_raw))
if (med_raw > 20) {
  cat("  -> umol/L detected; converting to mg/dL (divide by 88.42)\n")
  cr[, value := value / 88.42]
} else {
  cat("  -> already mg/dL; no conversion applied\n")
}
cr <- cr[value > 0.05 & value < 50]

ser <- merge(cr, episodes[, .(episode_id, patient_id, t0)],
             by = "patient_id", allow.cartesian = TRUE)
ser[, t_days := as.numeric(difftime(charttime, t0, units = "days"))]
ser <- ser[t_days >= -LOOKBACK_DAYS & t_days <= FOLLOWUP_DAYS]

series <- ser[, .(source = "HiRID", episode_id, patient_id, charttime,
                  t_days = round(t_days, 5), value = round(value, 4), unit = "mg/dL")]
setorder(series, episode_id, t_days)
fwrite(series, paste0(out_dir, "HiRID_creatinine_series.csv"))
cat(sprintf("  series written: %s rows\n", format(nrow(series), big.mark = ",")))

# ---------------------------------------------------------------- 3. checks
cat("\n=== Sanity checks (paste this into your reply) ===\n")
cat(sprintf("  median creatinine: %.2f mg/dL  (expect ~0.7-1.2 after conversion)\n", median(series$value)))
cat(sprintf("  episodes with >=1 result: %s of %s (%.1f%%)\n",
    format(uniqueN(series$episode_id), big.mark = ","),
    format(nrow(episodes), big.mark = ","),
    100*uniqueN(series$episode_id)/nrow(episodes)))
n <- series[, .N, by = episode_id]
cat(sprintf("  results per episode: median %d (IQR %d-%d), max %d\n",
    as.integer(median(n$N)), as.integer(quantile(n$N,.25)),
    as.integer(quantile(n$N,.75)), max(n$N)))
cat(sprintf("  timing: pre-t0 %.1f%% | 0-7d after %.1f%% | 7-90d after %.1f%%\n",
    100*mean(series$t_days < 0),
    100*mean(series$t_days > 0 & series$t_days <= 7),
    100*mean(series$t_days > 7 & series$t_days <= 90)))
cat("  NOTE: pre-t0 and day 7-90 percentages are expected to be near zero.\n")
cat("        HiRID has no pre-admission history and ends at ICU discharge.\n")
cat(sprintf("  covariates: age %.1f%% | sex %.1f%% | discharge status %.1f%%\n",
    100*mean(!is.na(episodes$age)), 100*mean(!is.na(episodes$sex)),
    100*mean(!is.na(episodes$discharge_status))))
cat("  discharge time and death time are NOT AVAILABLE in HiRID (columns left empty).\n")
cat("\nDone. Upload BOTH files. Do not open them in Excel.\n")
