rm(list = ls())

# AmsterdamUMCdb: complete creatinine series for every ICU/MC episode
# Output 1: creatinine_series.csv (one row per available creatinine result)
# Output 2: episode_covariates.csv (one row per admission, including episodes
#           with no creatinine measurement)
#
# Important source limitations:
# - AmsterdamUMCdb stores de-identified times as millisecond offsets rather than
#   real calendar datetimes. Therefore t0/timestamp/discharge_time/death_time are
#   exported as millisecond offsets on the database's common patient timeline.
# - Exact age is unavailable; age is the midpoint of age_group (80+ -> 85), and
#   age_group is retained.
# - ASA and hospital-level LOS are not available in the supplied AmsterdamUMCdb
#   tables. Required columns are retained as NA; ICU/MC LOS is available.
# - All creatinine results available anywhere in numericitems for the same
#   patient are considered, including results attached to another ICU admission
#   and pre-admission laboratory results. Pure outpatient results not present in
#   AmsterdamUMCdb cannot be recovered.

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------- parameters
# Every machine-specific setting can be overridden with an environment variable.
# Example (PowerShell):
# $env:AMDS_OUT_DIR = "path/to/amsterdam_output"
# $env:AMDS_HOST = "localhost"; $env:AMDS_PORT = "5432"
# $env:AMDS_DB = "AMDS"; $env:AMDS_USER = "postgres"
# $env:PSQL_PATH = "C:/Program Files/PostgreSQL/16/bin/psql.exe"
# Set PGPASSWORD outside this script if the server requires a password.
env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

local_psql <- "C:/AmsterdamUMCdb/bin/psql.exe"
PSQL <- env_or("PSQL_PATH", Sys.which("psql"))
if (!nzchar(PSQL) && file.exists(local_psql)) PSQL <- local_psql

HOST   <- env_or("AMDS_HOST", "localhost")
PORT   <- env_or("AMDS_PORT", "5433")
USER   <- env_or("AMDS_USER", "postgres")
DB     <- env_or("AMDS_DB", "AMDS")
SCHEMA <- env_or("AMDS_SCHEMA", "public")

OUT_DIR <- env_or("AMDS_OUT_DIR",
                  "path/to/amsterdam_output")
LOOKBACK_DAYS <- as.integer(env_or("LOOKBACK_DAYS", "365"))
FOLLOWUP_DAYS <- as.integer(env_or("FOLLOWUP_DAYS", "90"))
MS_PER_DAY <- 1000 * 60 * 60 * 24

# Serum/blood creatinine in AmsterdamUMCdb. Both are recorded as micromol/L.
CREAT_ITEMIDS <- c(6836L, 9941L)
UMOL_L_PER_MG_DL <- 88.4

if (!nzchar(PSQL) || !file.exists(PSQL)) {
  stop("psql not found. Put psql on PATH or set PSQL_PATH.")
}
if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", SCHEMA)) stop("Invalid AMDS_SCHEMA")
if (is.na(LOOKBACK_DAYS) || LOOKBACK_DAYS < 0) stop("Invalid LOOKBACK_DAYS")
if (is.na(FOLLOWUP_DAYS) || FOLLOWUP_DAYS < 0) stop("Invalid FOLLOWUP_DAYS")

admissions_table <- paste0(SCHEMA, ".admissions")
numericitems_table <- paste0(SCHEMA, ".numericitems")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

series_path <- file.path(OUT_DIR, "creatinine_series.csv")
episode_path <- file.path(OUT_DIR, "episode_covariates.csv")
log_path <- file.path(OUT_DIR, "extraction_log.txt")
manifest_path <- file.path(OUT_DIR, "reproducibility_manifest.txt")

sql_path <- tempfile(fileext = ".sql")
on.exit(unlink(sql_path), add = TRUE)

sql_quote_path <- function(x) gsub("'", "''", gsub("\\\\", "/", x))

run_psql_file <- function(path) {
  out <- system2(PSQL,
                 c("-h", HOST, "-p", PORT, "-U", USER, "-d", DB,
                   "-v", "ON_ERROR_STOP=1", "-f", shQuote(path)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (length(out)) write(out, file = log_path, append = TRUE)
  as.integer(status)
}

check_connection <- function() {
  out <- suppressWarnings(system2(
    PSQL,
    c("-h", HOST, "-p", PORT, "-U", USER, "-d", DB,
      "-v", "ON_ERROR_STOP=1", "-At", "-c", shQuote("SELECT 1;")),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L || !any(trimws(out) == "1")) {
    stop("Cannot connect to PostgreSQL or database is not ready. ",
         "Check AMDS_HOST/AMDS_PORT/AMDS_DB/AMDS_USER and PGPASSWORD.\n",
         paste(out, collapse = "\n"))
  }
}

check_connection()

cat("AmsterdamUMCdb creatinine export\n", file = log_path)
cat(sprintf("Started: %s\n", format(Sys.time(), tz = "Asia/Shanghai")),
    file = log_path, append = TRUE)
cat(sprintf("Window: -%d to +%d days around each ICU/MC admission t0\n",
            LOOKBACK_DAYS, FOLLOWUP_DAYS), file = log_path, append = TRUE)
cat(sprintf("Creatinine itemids: %s\n", paste(CREAT_ITEMIDS, collapse = ", ")),
    file = log_path, append = TRUE)
cat("Unit conversion: micromol/L divided by 88.4 = mg/dL\n",
    file = log_path, append = TRUE)
cat(sprintf("Database: host=%s port=%s db=%s schema=%s user=%s\n",
            HOST, PORT, DB, SCHEMA, USER), file = log_path, append = TRUE)

# ---------------------------------------------------------------- 1. episodes
episode_sql <- sprintf(
"COPY (
  SELECT
    'AmsterdamUMCdb'::text AS source,
    a.admissionid::text AS episode_id,
    a.patientid::text AS patient_id,
    a.admittedat::bigint AS t0,
    CASE a.agegroup
      WHEN '18-39' THEN 28.5
      WHEN '40-49' THEN 44.5
      WHEN '50-59' THEN 54.5
      WHEN '60-69' THEN 64.5
      WHEN '70-79' THEN 74.5
      WHEN '80+'   THEN 85.0
      ELSE NULL
    END::double precision AS age,
    CASE
      WHEN lower(a.gender) IN ('man', 'male') THEN 'M'
      WHEN lower(a.gender) IN ('vrouw', 'female') THEN 'F'
      ELSE a.gender
    END::text AS sex,
    NULL::double precision AS asa,
    a.urgency::integer AS emergency,
    NULL::double precision AS hosp_los_days,
    ((a.dischargedat - a.admittedat)::double precision / %s)::double precision
      AS icu_los_days,
    a.dischargedat::bigint AS discharge_time,
    a.dateofdeath::bigint AS death_time,
    ((a.dischargedat - a.admittedat)::double precision / %s)::double precision
      AS discharge_t_days,
    CASE WHEN a.dateofdeath IS NULL THEN NULL ELSE
      ((a.dateofdeath - a.admittedat)::double precision / %s)::double precision
    END AS death_t_days,
    a.agegroup::text AS age_group,
    a.location::text AS icu_location,
    a.specialty::text AS specialty,
    a.admissioncount::integer AS admission_count,
    a.admissionyeargroup::text AS admission_year_group,
    CASE
      WHEN a.dateofdeath IS NOT NULL
       AND a.dateofdeath >= a.admittedat
       AND a.dateofdeath <= a.dischargedat THEN 1 ELSE 0
    END::integer AS died_during_episode,
    'age is age-group midpoint; ASA and hospital LOS unavailable'::text
      AS source_specific_note
  FROM %s a
  ORDER BY a.admissionid
) TO '%s' WITH (FORMAT CSV, HEADER TRUE, ENCODING 'UTF8');",
as.character(MS_PER_DAY), as.character(MS_PER_DAY), as.character(MS_PER_DAY),
admissions_table,
sql_quote_path(episode_path))

writeLines(c("\\set ON_ERROR_STOP on", episode_sql), sql_path, useBytes = TRUE)
cat("Exporting episode_covariates.csv ...\n")
status <- run_psql_file(sql_path)
if (!identical(status, 0L)) stop("Episode export failed; see ", log_path)

# ---------------------------------------------------------------- 2. series
# measuredat and admittedat use the same de-identified millisecond timeline.
# Joining by patient_id (not admission_id) captures available measurements from
# pre-admission periods and other/repeated admissions within each episode window.
series_sql <- sprintf(
"COPY (
  WITH episodes AS (
    SELECT
      a.admissionid,
      a.patientid,
      a.admittedat AS t0
    FROM %s a
  ),
  all_patient_creatinine AS (
    SELECT
      src.patientid,
      n.admissionid AS source_admissionid,
      n.itemid,
      n.measuredat,
      n.value AS creatinine_umol_l,
      n.unit AS original_unit
    FROM %s n
    INNER JOIN %s src
      ON src.admissionid = n.admissionid
    WHERE n.itemid IN (%s)
      AND n.value IS NOT NULL
      AND n.measuredat IS NOT NULL
      -- Match the supplied templates: retain every positive result below
      -- 50 mg/dL. In native units this is 50 * 88.4 = 4420 micromol/L.
      AND n.value > 0
      AND n.value < (50 * %s)
  )
  SELECT
    'AmsterdamUMCdb'::text AS source,
    ep.admissionid::text AS episode_id,
    ep.patientid::text AS patient_id,
    cr.measuredat::bigint AS timestamp,
    round(((cr.measuredat - ep.t0)::numeric / %s), 5) AS t_days,
    (cr.creatinine_umol_l / %s)::double precision AS value,
    'mg/dL'::text AS unit,
    cr.source_admissionid::text AS source_admission_id,
    cr.itemid::bigint AS source_itemid,
    cr.creatinine_umol_l::double precision AS original_value,
    cr.original_unit::text AS original_unit
  FROM episodes ep
  INNER JOIN all_patient_creatinine cr
    ON cr.patientid = ep.patientid
   AND cr.measuredat >= ep.t0 - (%s::bigint * %s::bigint)
   AND cr.measuredat <= ep.t0 + (%s::bigint * %s::bigint)
  ORDER BY ep.admissionid, cr.measuredat, cr.source_admissionid, cr.itemid
) TO '%s' WITH (FORMAT CSV, HEADER TRUE, ENCODING 'UTF8');",
admissions_table, numericitems_table, admissions_table,
paste(CREAT_ITEMIDS, collapse = ", "), as.character(UMOL_L_PER_MG_DL),
as.character(MS_PER_DAY), as.character(UMOL_L_PER_MG_DL),
as.character(LOOKBACK_DAYS), as.character(MS_PER_DAY),
as.character(FOLLOWUP_DAYS), as.character(MS_PER_DAY),
sql_quote_path(series_path))

writeLines(c("\\set ON_ERROR_STOP on", series_sql), sql_path, useBytes = TRUE)
cat("Exporting creatinine_series.csv ...\n")
status <- run_psql_file(sql_path)
if (!identical(status, 0L)) stop("Series export failed; see ", log_path)

# ---------------------------------------------------------------- 3. checks
cat("Running sanity checks ...\n")
episodes <- fread(episode_path, colClasses = list(character = c("episode_id", "patient_id")))
series <- fread(series_path,
                colClasses = list(character = c("episode_id", "patient_id",
                                                "source_admission_id")))

required_episode <- c("episode_id", "patient_id", "t0", "age", "sex", "asa",
                      "emergency", "hosp_los_days", "icu_los_days",
                      "discharge_time", "death_time", "discharge_t_days",
                      "death_t_days")
required_series <- c("episode_id", "patient_id", "timestamp", "t_days", "value", "unit")
stopifnot(all(required_episode %in% names(episodes)))
stopifnot(all(required_series %in% names(series)))
stopifnot(nrow(episodes) == uniqueN(episodes$episode_id))
stopifnot(all(series$t_days >= -LOOKBACK_DAYS & series$t_days <= FOLLOWUP_DAYS))
stopifnot(all(series$unit == "mg/dL"))

n_with <- uniqueN(series$episode_id)
n_without <- nrow(episodes) - n_with
n_per <- series[, .N, by = episode_id]

qc <- c(
  sprintf("Completed: %s", format(Sys.time(), tz = "Asia/Shanghai")),
  sprintf("Episodes: %s", format(nrow(episodes), big.mark = ",")),
  sprintf("Unique patients: %s", format(uniqueN(episodes$patient_id), big.mark = ",")),
  sprintf("Creatinine rows: %s", format(nrow(series), big.mark = ",")),
  sprintf("Episodes with >=1 result: %s (%.2f%%)", format(n_with, big.mark = ","),
          100 * n_with / nrow(episodes)),
  sprintf("Episodes with zero results (retained in covariate file): %s",
          format(n_without, big.mark = ",")),
  sprintf("Median creatinine: %.3f mg/dL", median(series$value, na.rm = TRUE)),
  sprintf("Results/episode among measured: median %d, IQR %d-%d, max %d",
          as.integer(median(n_per$N)), as.integer(quantile(n_per$N, .25)),
          as.integer(quantile(n_per$N, .75)), max(n_per$N)),
  sprintf("Pre-t0 rows: %s (%.2f%%)", format(sum(series$t_days < 0), big.mark = ","),
          100 * mean(series$t_days < 0)),
  sprintf("Day 7-90 rows: %s (%.2f%%)",
          format(sum(series$t_days > 7 & series$t_days <= 90), big.mark = ","),
          100 * mean(series$t_days > 7 & series$t_days <= 90)),
  "Age is approximated from age-group midpoint (80+ -> 85).",
  "ASA and hospital-level LOS are unavailable and exported as NA.",
  "Death date is day-level/de-identified; death_t_days may be slightly negative when ICU admission occurred later on the recorded death date.",
  "All available cross-admission/preadmission creatinine records were searched.",
  "Pure outpatient data absent from AmsterdamUMCdb cannot be exported."
)
write(qc, file = log_path, append = TRUE)

manifest <- c(
  "AmsterdamUMCdb creatinine export reproducibility manifest",
  paste0("generated_at=", format(Sys.time(), tz = "Asia/Shanghai")),
  paste0("r_version=", R.version.string),
  paste0("data_table_version=", as.character(packageVersion("data.table"))),
  paste0("psql_version=", system2(PSQL, "--version", stdout = TRUE)[1]),
  paste0("database_host=", HOST),
  paste0("database_port=", PORT),
  paste0("database_name=", DB),
  paste0("database_schema=", SCHEMA),
  paste0("episode_definition=one ICU/MC admission; t0=admissions.admittedat"),
  paste0("lookback_days=", LOOKBACK_DAYS),
  paste0("followup_days=", FOLLOWUP_DAYS),
  paste0("creatinine_itemids=", paste(CREAT_ITEMIDS, collapse = ",")),
  paste0("conversion=original_micromol_per_litre/", UMOL_L_PER_MG_DL),
  paste0("episode_rows=", nrow(episodes)),
  paste0("unique_patients=", uniqueN(episodes$patient_id)),
  paste0("series_rows=", nrow(series)),
  paste0("episode_covariates_md5=", unname(tools::md5sum(episode_path))),
  paste0("creatinine_series_md5=", unname(tools::md5sum(series_path)))
)
writeLines(manifest, manifest_path, useBytes = TRUE)
cat(paste(qc, collapse = "\n"), "\n")
cat("\nCreated:\n", episode_path, "\n", series_path, "\n", sep = "")
