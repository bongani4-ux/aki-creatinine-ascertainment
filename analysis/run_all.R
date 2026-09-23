# =============================================================================
# run_all.R  -  run the full analysis. Each source runs in its own R process to
# bound memory, then results are combined.
#
# Usage:  Rscript run_all.R                 all sources
#         Rscript run_all.R MIMIC-IV PIC    named sources only
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source("config.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2 && args[1] == "--one") {
  source("pipeline.R")
  k <- args[2]; row <- SOURCES[SOURCES$key == k, ]
  p <- source_paths(k)
  res <- analyse_source(k, p[["series"]], p[["covar"]], row$setting, row$clock, row$offset_unit)
  dir.create(file.path(RESULTS_DIR, "per_source"), recursive = TRUE, showWarnings = FALSE)
  write_json(res, file.path(RESULTS_DIR, "per_source", paste0(k, ".json")),
             auto_unbox = TRUE, digits = NA, null = "null")
  cat(sprintf("%-15s n=%9s  rule=%s\n", k, format(res$measure$n, big.mark = ","),
              res$measure$first_episode_rule))
  quit(save = "no")
}

keys <- if (length(args)) args else SOURCES$key
rscript <- file.path(R.home("bin"), "Rscript")
for (k in keys) {
  status <- system2(rscript, c("run_all.R", "--one", shQuote(k)))
  if (status != 0) stop("analysis failed for ", k)
}

# ---- combine
per <- file.path(RESULTS_DIR, "per_source")
files <- file.path(per, paste0(SOURCES$key, ".json")); files <- files[file.exists(files)]
res <- lapply(files, read_json, simplifyVector = TRUE)
tbl <- function(part) rbindlist(lapply(res, function(x) {
  v <- x[[part]]; if (is.null(v) || !length(v)) return(NULL)
  as.data.table(lapply(v, function(e) if (is.null(e)) NA else e)) }), fill = TRUE)
fwrite(tbl("measure"),  file.path(RESULTS_DIR, "measurement_and_provenance.csv"))
fwrite(tbl("rules"),    file.path(RESULTS_DIR, "ascertainment_by_baseline_rule.csv"))
fwrite(tbl("misclass"), file.path(RESULTS_DIR, "baseline_misclassification.csv"))
timing <- setNames(lapply(res, `[[`, "timing"), sapply(res, function(x) x$measure$db))
write_json(timing, file.path(RESULTS_DIR, "timing.json"), auto_unbox = TRUE, digits = NA, null = "null")
cat(sprintf("combined %d sources -> %s/\n", length(res), RESULTS_DIR))
