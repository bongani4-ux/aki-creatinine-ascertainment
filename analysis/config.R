# =============================================================================
# config.R  -  source registry. Point DATA_DIR at the folder holding the exports.
# Each source needs <KEY>_creatinine_series.csv and <KEY>_episode_covariates.csv.
# CLOCK names the column holding the clock time of each draw, where it exists.
# OFFSET_UNIT gives the unit of numeric discharge/death offsets relative to t0
# (eICU-CRD minutes, SICdb and VitalDB seconds, AmsterdamUMCdb milliseconds).
# =============================================================================
DATA_DIR    <- "data"
RESULTS_DIR <- "results"

# Remove exact duplicate results (same episode, same time, same value) before
# counting. They are repeated records of one draw, not separate measurements
# (60% of HiRID rows, 8.5% of INSPIRE, ~7% of SICdb and PIC in the exports).
DEDUPLICATE <- TRUE

SOURCES <- data.frame(
  key = c("MIMIC-IV", "eICU-CRD", "AmsterdamUMCdb", "SICdb", "NWICU", "HiRID", "Zigong",
          "PIC", "VitalDB", "INSPIRE", "MOVER", "NJF", "ZD", "AHMU"),
  setting = c(rep("Critical care", 7), "Critical care (paed)",
              rep("Perioperative", 3), rep("Surgical (inst)", 3)),
  clock = c("charttime", NA, NA, NA, "charttime", "charttime", NA,
            "charttime", NA, NA, "timestamp", NA, NA, NA),
  # unit of numeric discharge_time / death_time offsets, where not datetimes
  offset_unit = c(NA, "minutes", "milliseconds", "seconds", NA, NA, NA,
                  NA, "seconds", NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE)

source_paths <- function(key)
  c(series = file.path(DATA_DIR, paste0(key, "_creatinine_series.csv")),
    covar  = file.path(DATA_DIR, paste0(key, "_episode_covariates.csv")))
