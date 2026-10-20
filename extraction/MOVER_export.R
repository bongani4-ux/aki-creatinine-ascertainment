rm(list = ls())

###############################################################################
# MOVER (EPIC subset) - longitudinal creatinine export
#
# Input files (MOVER EPIC_EMR/EMR/):
#   1. patient_information.csv
#   2. patient_labs.csv
#
# Output files (written to out_dir):
#   1. MOVER_episode_covariates.csv
#   2. MOVER_creatinine_series.csv
#
# patient_labs.csv (~3.5 GB) is read in a single fread() call, selecting 6 columns.
###############################################################################

library(data.table)


###############################################################################
# 1. Paths and parameters
###############################################################################

# Folder holding patient_information.csv and patient_labs.csv, and output folder.
script_dir <- "path/to/mover/EPIC_EMR/EMR"     # input folder
out_dir    <- "path/to/export_series/"         # output folder
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

info_file <-
  file.path(
    script_dir,
    "patient_information.csv"
  )

labs_file <-
  file.path(
    script_dir,
    "patient_labs.csv"
  )

# Export window: 365 days before to 90 days after t0, inclusive.
LOOKBACK_DAYS <- 365
FOLLOWUP_DAYS <- 90

# Serum/plasma (2160-0) and blood (38483-4) creatinine LOINC codes.
# A warning is raised if Lab Name = Creatinine appears with any other code.
CREAT_CODES <- c(
  "2160-0",
  "38483-4"
)

# All arithmetic in UTC so daylight-saving changes do not distort time differences.
TIME_ZONE <- "UTC"


###############################################################################
# 2. Helpers
###############################################################################

# Parse the date/datetime formats that occur in MOVER (m/d/yy H:MM in practice).
parse_mover_datetime <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x[x == ""] <- NA_character_

  mover_formats <- c(
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%m/%d/%y %H:%M",
    "%m/%d/%Y %H:%M",
    "%Y/%m/%d %H:%M:%S",
    "%Y-%m-%d",
    "%m/%d/%y",
    "%m/%d/%Y"
  )

  parsed_time <- as.POSIXct(
    rep(
      NA_real_,
      length(x)
    ),
    origin = "1970-01-01",
    tz = TIME_ZONE
  )

  for (current_format in mover_formats) {

    rows_to_parse <- which(
      is.na(parsed_time) &
        !is.na(x)
    )

    if (!length(rows_to_parse)) {
      break
    }

    current_result <- suppressWarnings(
      as.POSIXct(
        x[rows_to_parse],
        format = current_format,
        tz = TIME_ZONE
      )
    )

    successfully_parsed <-
      !is.na(current_result)

    if (any(successfully_parsed)) {

      parsed_time[
        rows_to_parse[successfully_parsed]
      ] <-
        current_result[successfully_parsed]
    }
  }

  parsed_time
}


# Write times as YYYY-MM-DDTHH:MM:SS.
iso_time <- function(x) {

  ifelse(
    is.na(x),
    NA_character_,
    format(
      x,
      "%Y-%m-%dT%H:%M:%S",
      tz = TIME_ZONE
    )
  )
}


# Stop if a required column is missing.
require_columns <- function(
    column_names,
    required_columns,
    table_name) {

  missing_columns <- setdiff(
    required_columns,
    column_names
  )

  if (length(missing_columns)) {

    stop(
      table_name,
      " is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
}


# Normalise how the creatinine unit is written.
normalize_unit <- function(x) {

  x <- enc2utf8(
    ifelse(
      is.na(x),
      "",
      x
    )
  )

# Map the different micro signs (µ, μ, Μ) to U.
  x <- gsub(
    "[\u00B5\u03BC\u039C]",
    "U",
    x
  )

# Remove whitespace and upper-case.
  toupper(
    gsub(
      "[[:space:]]+",
      "",
      x
    )
  )
}


###############################################################################
# 3. Episodes (one row per operation)
###############################################################################

cat(
  "=== Building episode table ===\n"
)

info <- fread(
  info_file,
  colClasses = "character",
  check.names = FALSE,
  showProgress = TRUE
)

raw_episode_rows <- nrow(info)

required_info_columns <- c(
  "LOG_ID",
  "MRN",
  "DISCH_DISP",
  "HOSP_ADMSN_TIME",
  "HOSP_DISCH_TIME",
  "ICU_ADMIN_FLAG",
  "SURGERY_DATE",
  "BIRTH_DATE",
  "SEX",
  "ASA_RATING_C",
  "PRIMARY_PROCEDURE_NM",
  "IN_OR_DTTM",
  "AN_START_DATETIME"
)

require_columns(
  names(info),
  required_info_columns,
  "patient_information.csv"
)

# Trim episode and patient identifiers.
info[, LOG_ID := trimws(LOG_ID)]
info[, MRN := trimws(MRN)]

if (
  anyNA(info$LOG_ID) ||
    any(info$LOG_ID == "")
) {

  stop(
    "LOG_ID is missing for at least one source row."
  )
}

if (
  anyNA(info$MRN) ||
    any(info$MRN == "")
) {

  stop(
    "MRN is missing for at least one source row."
  )
}

# Drop rows that are exact duplicates across all fields.
info <- unique(info)

# LOG_ID recurs on a few different surgery dates, so an episode is
# identified by LOG_ID + SURGERY_DATE.
patient_conflicts <- info[
  ,
  .(
    n_patients__ = uniqueN(MRN)
  ),
  by = .(
    LOG_ID,
    SURGERY_DATE
  )
][
  n_patients__ > 1L
]

if (nrow(patient_conflicts)) {

  warning(
    nrow(patient_conflicts),
    " source episode(s) map to more than one MRN representation. ",
    "No MRN value was changed; the first source row will be retained."
  )
}

# Keep source row order.
info[, source_order__ := .I]

setorder(
  info,
  LOG_ID,
  SURGERY_DATE,
  source_order__
)

# Keep the first source row per LOG_ID + SURGERY_DATE.
info <- info[
  ,
  .SD[1L],
  by = .(
    LOG_ID,
    SURGERY_DATE
  )
]

# Parse dates and datetimes.
info[, `:=`(
  surgery_dt__ =
    parse_mover_datetime(
      SURGERY_DATE
    ),

  in_or_dt__ =
    parse_mover_datetime(
      IN_OR_DTTM
    ),

  an_start_dt__ =
    parse_mover_datetime(
      AN_START_DATETIME
    ),

  admission_dt__ =
    parse_mover_datetime(
      HOSP_ADMSN_TIME
    ),

  discharge_dt__ =
    parse_mover_datetime(
      HOSP_DISCH_TIME
    )
)]

# t0 = entry to the operating room; if missing, anaesthesia start; if that is also
# missing, SURGERY_DATE, which is midnight of the operation day (see t0_source).
info[, t0__ := in_or_dt__]
info[, t0_source__ := fifelse(!is.na(in_or_dt__), "in_or",
                       fifelse(!is.na(an_start_dt__), "anaesthesia_start", "surgery_date_midnight"))]

info[
  is.na(t0__),
  t0__ := an_start_dt__
]

info[
  is.na(t0__),
  t0__ := surgery_dt__
]

if (anyNA(info$t0__)) {

  stop(
    "t0 could not be derived for at least one episode."
  )
}

# Guard against a date format being misread (e.g. a 4-digit year read by %y).
t0_year__ <- as.integer(format(info$t0__, "%Y"))
if (any(t0_year__ < 2000 | t0_year__ > 2030)) {
  stop("Implausible t0 year(s): ", paste(sort(unique(t0_year__[t0_year__ < 2000 | t0_year__ > 2030])),
                                         collapse = ", "), ". Check parse_mover_datetime().")
}
cat("t0 source:\n"); print(table(info$t0_source__))

# Where one LOG_ID has several episodes, order by t0 and number them.
setorder(
  info,
  LOG_ID,
  t0__,
  MRN
)

info[, `:=`(
  episode_n__ = .N,
  episode_seq__ = seq_len(.N)
), by = LOG_ID]

info[, episode_id__ := fifelse(
  episode_n__ > 1L,
  paste0(
    LOG_ID,
    "__",
    episode_seq__
  ),
  LOG_ID
)]

# Hospital length of stay.
info[, hosp_los__ := as.numeric(
  difftime(
    discharge_dt__,
    admission_dt__,
    units = "days"
  )
)]

# Hospital LOS below 0 or above 365 days set to missing.
info[
  !between(
    hosp_los__,
    0,
    365
  ),
  hosp_los__ := NA_real_
]

# In MOVER EPIC, BIRTH_DATE holds the de-identified age in years.
info[, age__ := suppressWarnings(
  as.numeric(BIRTH_DATE)
)]

info[
  !between(
    age__,
    0,
    110
  ),
  age__ := NA_real_
]

# ASA class: keep 1-6.
info[, asa__ := suppressWarnings(
  as.integer(
    as.numeric(ASA_RATING_C)
  )
)]

info[
  !between(
    asa__,
    1L,
    6L
  ),
  asa__ := NA_integer_
]

# Emergency only where the procedure name explicitly contains 'emerg'.
procedure_text__ <- ifelse(
  is.na(info$PRIMARY_PROCEDURE_NM),
  "",
  info$PRIMARY_PROCEDURE_NM
)

procedure_text__ <- tolower(
  trimws(procedure_text__)
)

explicit_emergency__ <- grepl(
  "emerg",
  procedure_text__,
  fixed = TRUE
)

# ICU admission flag.
icu_flag__ <- ifelse(
  is.na(info$ICU_ADMIN_FLAG),
  "",
  info$ICU_ADMIN_FLAG
)

icu_flag__ <- tolower(
  trimws(icu_flag__)
)

# Final episode covariate table.
episodes <- info[, .(
  episode_id = episode_id__,
  patient_id = MRN,
  t0 = iso_time(t0__),
  t0_source = t0_source__,
  age = age__,

  sex = fcase(
    tolower(trimws(SEX)) == "male",
    "M",

    tolower(trimws(SEX)) == "female",
    "F",

    default = NA_character_
  ),

  ASA = asa__,

  emergency = fifelse(
    explicit_emergency__,
    1L,
    NA_integer_
  ),

  hosp_los_days = round(
    hosp_los__,
    5
  ),

  # 0 if never admitted to ICU; missing if admitted (no ICU in/out times in MOVER).
  icu_los_days = fifelse(
    icu_flag__ == "no",
    0,
    NA_real_
  ),

  discharge_time =
    iso_time(discharge_dt__),

  # Only in-hospital death is identifiable; discharge time is used as the death time.
  death_time = fifelse(
    tolower(trimws(DISCH_DISP)) == "expired",
    iso_time(discharge_dt__),
    NA_character_
  )
)]

if (
  uniqueN(episodes$episode_id) !=
    nrow(episodes)
) {

  stop(
    "episode_id is not unique after source de-duplication."
  )
}

# Write the episode table before reading labs, so patients with no creatinine
# are kept in it.
episode_output <- file.path(
  out_dir,
  "MOVER_episode_covariates.csv"
)

fwrite(
  episodes,
  episode_output,
  na = ""
)

# Lookup used to link creatinine results to index episodes.
episode_lookup <- info[, .(
  episode_id = episode_id__,
  patient_id = MRN,
  t0__
)]


###############################################################################
# 4. Laboratory data
###############################################################################

cat(
  "=== Reading patient_labs.csv in one fread call ===\n"
)

needed_lab_columns <- c(
  "MRN",
  "Lab Code",
  "Lab Name",
  "Observation Value",
  "Measurement Units",
  "Collection Datetime"
)

# Read only the 6 columns needed.
#
labs <- fread(
  labs_file,
  select = needed_lab_columns,
  colClasses = "character",
  check.names = FALSE,
  showProgress = TRUE
)

require_columns(
  names(labs),
  needed_lab_columns,
  "patient_labs.csv"
)

raw_lab_rows <- nrow(labs)

# Clean the fields used for filtering and joining.
labs[, MRN := trimws(MRN)]

labs[, lab_name__ := tolower(
  trimws(
    fifelse(
      is.na(`Lab Name`),
      "",
      `Lab Name`
    )
  )
)]

labs[, lab_code__ := trimws(
  fifelse(
    is.na(`Lab Code`),
    "",
    `Lab Code`
  )
)]

# Codes that actually occur with Lab Name = Creatinine.
observed_creatinine_codes <- sort(
  unique(
    labs[
      lab_name__ == "creatinine" &
        lab_code__ != "",
      lab_code__
    ]
  )
)

# Warn about other codes labelled Creatinine but do not include them
# automatically, so an unverified item cannot enter the data.
unselected_creatinine_codes <- setdiff(
  observed_creatinine_codes,
  CREAT_CODES
)

if (length(unselected_creatinine_codes)) {

  warning(
    "Other Lab Code value(s) were found for Lab Name = Creatinine but were not selected: ",
    paste(
      unselected_creatinine_codes,
      collapse = ", "
    ),
    ". Review these codes before final submission."
  )
}

named_creatinine_rows <- sum(
  labs$lab_name__ == "creatinine"
)

blood_code_rows <- sum(
  labs$lab_name__ == "creatinine" &
    labs$lab_code__ %chin% CREAT_CODES
)

# Select rows where:
#   1. Lab Name is Creatinine;
#   2. Lab Code is in CREAT_CODES;
#   3. MRN belongs to the episode table.
cr <- labs[
  lab_name__ == "creatinine" &
    lab_code__ %chin% CREAT_CODES &
    MRN %chin% unique(episode_lookup$patient_id),
  .(
    patient_id = MRN,
    value_text = `Observation Value`,
    unit_raw = `Measurement Units`,
    timestamp_text = `Collection Datetime`
  )
]

cohort_blood_rows <- nrow(cr)

# Parse the collection datetime.
cr[, timestamp :=
  parse_mover_datetime(
    timestamp_text
  )
]

# Result value to numeric.
cr[, value := suppressWarnings(
  as.numeric(value_text)
)]

# Normalise unit text.
cr[, unit_key__ :=
  normalize_unit(unit_raw)
]

# Flag values recorded in umol/L.
cr[, converted__ :=
  unit_key__ == "UMOL/L"
]

# Only mg/dL and umol/L are supported.
cr[, supported_unit__ :=
  unit_key__ %chin%
    c(
      "MG/DL",
      "UMOL/L"
    )
]

converted_input_rows <- sum(
  cr$converted__,
  na.rm = TRUE
)

# umol/L / 88.4 -> mg/dL.
cr[
  converted__ == TRUE,
  value := value / 88.4
]

# Drop rows with an invalid patient, time, value or unit,
# and keep 0 < value < 50 mg/dL.
cr <- cr[
  supported_unit__ &
    !is.na(patient_id) &
    patient_id != "" &
    !is.na(timestamp) &
    !is.na(value) &
    value > 0 &
    value < 50
]

valid_rows_before_deduplication <- nrow(cr)

if (!nrow(cr)) {

  stop(
    "No valid serum/plasma or whole-blood creatinine rows were found."
  )
}

# One row per patient, time, value and unit.
cr[, unit := "mg/dL"]

cr <- cr[
  ,
  .(
    converted__ = any(converted__)
  ),
  by = .(
    patient_id,
    timestamp,
    value,
    unit
  )
]

unique_lab_rows <- nrow(cr)


###############################################################################
# 5. Link every result to every episode of the same patient
###############################################################################

# Join on patient (MRN), not LOG_ID, so outpatient and readmission results
# within the window are included.
ser <- merge(
  cr,
  episode_lookup,
  by = "patient_id",
  allow.cartesian = TRUE
)

# Days from this episode's t0.
ser[, t_days := as.numeric(
  difftime(
    timestamp,
    t0__,
    units = "days"
  )
)]

# Keep -365 <= t_days <= 90, inclusive.
ser <- ser[
  between(
    t_days,
    -LOOKBACK_DAYS,
    FOLLOWUP_DAYS
  )
]

converted_export_rows <- sum(
  ser$converted__,
  na.rm = TRUE
)

# Final series columns.
series <- ser[, .(
  episode_id,
  patient_id,
  timestamp = iso_time(timestamp),
  t_days = round(
    t_days,
    5
  ),
  value = round(
    value,
    6
  ),
  unit
)]

series <- unique(series)

setorder(
  series,
  episode_id,
  t_days,
  timestamp,
  value
)

series_output <- file.path(
  out_dir,
  "MOVER_creatinine_series.csv"
)

fwrite(
  series,
  series_output,
  na = ""
)


###############################################################################
# 6. Validation
###############################################################################

expected_episode_columns <- c(
  "episode_id",
  "patient_id",
  "t0",
  "t0_source",
  "age",
  "sex",
  "ASA",
  "emergency",
  "hosp_los_days",
  "icu_los_days",
  "discharge_time",
  "death_time"
)

expected_series_columns <- c(
  "episode_id",
  "patient_id",
  "timestamp",
  "t_days",
  "value",
  "unit"
)

if (
  !identical(
    names(episodes),
    expected_episode_columns
  )
) {

  stop(
    "MOVER_episode_covariates.csv columns do not match the required schema."
  )
}

if (
  !identical(
    names(series),
    expected_series_columns
  )
) {

  stop(
    "MOVER_creatinine_series.csv columns do not match the required schema."
  )
}

if (
  nrow(series) &&
    !all(
      between(
        series$t_days,
        -LOOKBACK_DAYS,
        FOLLOWUP_DAYS
      )
    )
) {

  stop(
    "creatinine_series.csv contains rows outside the requested time window."
  )
}

if (
  nrow(series) &&
    !all(series$unit == "mg/dL")
) {

  stop(
    "creatinine_series.csv contains a non-mg/dL unit."
  )
}

if (
  !all(
    series$episode_id %chin%
      episodes$episode_id
  )
) {

  stop(
    "creatinine_series.csv contains an unknown episode_id."
  )
}

episodes_with_result <- uniqueN(
  series$episode_id
)

episodes_without_result <-
  nrow(episodes) -
  episodes_with_result

median_creatinine <- if (nrow(series)) {
  median(series$value)
} else {
  NA_real_
}


###############################################################################
# 7. Summary
###############################################################################

cat(
  "\n=== Sanity checks ===\n"
)

cat(sprintf(
  "raw episode rows: %s\n",
  format(
    raw_episode_rows,
    big.mark = ","
  )
))

cat(sprintf(
  "unique episodes: %s; patients: %s\n",
  format(
    nrow(episodes),
    big.mark = ","
  ),
  format(
    uniqueN(episodes$patient_id),
    big.mark = ","
  )
))

cat(sprintf(
  "laboratory rows read: %s\n",
  format(
    raw_lab_rows,
    big.mark = ","
  )
))

cat(sprintf(
  "Lab Name = Creatinine rows: %s; selected-code rows: %s\n",
  format(
    named_creatinine_rows,
    big.mark = ","
  ),
  format(
    blood_code_rows,
    big.mark = ","
  )
))

cat(sprintf(
  "selected cohort rows before validation: %s\n",
  format(
    cohort_blood_rows,
    big.mark = ","
  )
))

cat(sprintf(
  "valid cohort lab rows before/after de-duplication: %s / %s\n",
  format(
    valid_rows_before_deduplication,
    big.mark = ","
  ),
  format(
    unique_lab_rows,
    big.mark = ","
  )
))

cat(sprintf(
  "series rows: %s\n",
  format(
    nrow(series),
    big.mark = ","
  )
))

cat(sprintf(
  "episodes with >=1 result: %s; without a result: %s\n",
  format(
    episodes_with_result,
    big.mark = ","
  ),
  format(
    episodes_without_result,
    big.mark = ","
  )
))

cat(sprintf(
  "median creatinine: %.2f mg/dL\n",
  median_creatinine
))

cat(sprintf(
  paste0(
    "unit conversion: %s input rows and %s exported ",
    "episode-result rows converted from umol/L to mg/dL\n"
  ),
  format(
    converted_input_rows,
    big.mark = ","
  ),
  format(
    converted_export_rows,
    big.mark = ","
  )
))

cat(sprintf(
  "time-window check: min %.5f days; max %.5f days\n",
  if (nrow(series)) {
    min(series$t_days)
  } else {
    NA_real_
  },
  if (nrow(series)) {
    max(series$t_days)
  } else {
    NA_real_
  }
))

cat(
  "\nDone. Keep and send BOTH files separately:\n"
)

cat(
  "  ",
  episode_output,
  "\n"
)

cat(
  "  ",
  series_output,
  "\n"
)
