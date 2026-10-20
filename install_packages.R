# Install the R packages the analysis needs. Run once:  Rscript install_packages.R
# Versions the paper was produced with are listed in README.md.
pkgs <- c("data.table", "jsonlite", "ggplot2", "patchwork", "scales", "ragg")
# The extraction scripts additionally use duckdb + DBI (MIMIC-IV, NWICU) and,
# for HiRID parquet files, arrow. Install those only if you run the exports.
extraction <- c("duckdb", "DBI")
need <- setdiff(c(pkgs, extraction), rownames(installed.packages()))
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
for (p in pkgs) cat(sprintf("%-11s %s\n", p, as.character(packageVersion(p))))
