rm(list = ls())

# NJF creatinine export. The Python helper performs one robust streaming scan
# of the large CSV and writes the two requested result files.
HOSPITAL <- "NJF"
LOOKBACK_DAYS <- 365L
FOLLOWUP_DAYS <- 90L

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (is.na(script_file)) stop("Run this file with Rscript so its directory can be resolved.")
script_dir <- dirname(script_file)

python <- "path/to/python.exe"
duckdb_dir <- "path/to/duckdb_temp"
input_root <- "path/to/institutional_input"
output_root <- "path/to/institutional_output"
helper <- file.path(script_dir, "export_creatinine.py")
qa_output <- file.path(script_dir, "..", "qa", paste0(HOSPITAL, "_qa_summary.json"))

args <- c(
  helper,
  "--hospital", HOSPITAL,
  "--input-root", input_root,
  "--output-root", output_root,
  "--duckdb-dir", duckdb_dir,
  "--qa-output", qa_output,
  "--lookback-days", LOOKBACK_DAYS,
  "--followup-days", FOLLOWUP_DAYS
)
output <- system2(python, args = shQuote(args), stdout = TRUE, stderr = TRUE)
status <- attr(output, "status")
cat(paste(output, collapse = "\n"), "\n")
if (!is.null(status) && status != 0L) stop("NJF export failed with status ", status)
