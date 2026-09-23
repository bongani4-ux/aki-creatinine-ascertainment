# =============================================================================
# tables.R  -  build the numeric manuscript tables (CSV) from results/.
# Usage:  Rscript tables.R        writes tables/*.csv
# Table 1 (literature) is a curated evidence table and is not generated here.
# Numbering follows the manuscript: Table 2 measurement intensity, Table 3
# baseline provenance, Table 4 baseline misclassification, Table 5 ascertainment
# by baseline rule; Supplementary Tables S1 timing, S2-S3 undefined vs unmeasured.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source("config.R")
OUT <- "tables"; dir.create(OUT, showWarnings = FALSE)
ORDER <- SOURCES$key
COUNTRY <- c("MIMIC-IV" = "USA", "eICU-CRD" = "USA", "AmsterdamUMCdb" = "Netherlands",
             "SICdb" = "Austria", "NWICU" = "USA", "HiRID" = "Switzerland", "Zigong" = "China",
             "PIC" = "China", "VitalDB" = "South Korea", "INSPIRE" = "South Korea", "MOVER" = "USA",
             "NJF" = "China", "ZD" = "China", "AHMU" = "China")
f1 <- function(x) sprintf("%.1f", x)
fn <- function(x) format(round(x), big.mark = ",", trim = TRUE)
ord <- function(d) d[order(match(db, ORDER))]

R  <- ord(fread(file.path(RESULTS_DIR, "measurement_and_provenance.csv")))
C3 <- ord(fread(file.path(RESULTS_DIR, "ascertainment_by_baseline_rule.csv")))
M  <- ord(fread(file.path(RESULTS_DIR, "baseline_misclassification.csv")))
TM <- read_json(file.path(RESULTS_DIR, "timing.json"), simplifyVector = TRUE)
out <- function(d, f) fwrite(d, file.path(OUT, f))

out(R[, .(`Data source` = db, Setting = setting, Country = COUNTRY[db], `Patients, n` = fn(n),
          `Creatinine results, n` = fn(results),
          `Results per patient, median (IQR)` = sprintf("%d (%d–%d)", med, q1, q3),
          `95th centile` = fn(p95), Maximum = fn(mx), `Never measured, %` = f1(zero))],
    "Table2_measurement_intensity.csv")

out(R[, .(`Data source` = db, `Tier 1 (7–30 d), %` = f1(T1), `Tier 2 (31–90 d), %` = f1(T2),
          `Tier 3 (91–180 d), %` = f1(T3), `Tier 4 (181–365 d), %` = f1(T4),
          `Hierarchy-compliant, %` = f1(hier), `Tier 0 proximal (<7 d), %` = f1(t0),
          `No pre-index value, %` = f1(nopre))], "Table3_baseline_provenance.csv")

out(M[, .(`Data source` = db, `Patients, n` = fn(n),
          `Hierarchy baseline, median mg/dL (µmol/L)` = sprintf("%.2f (%.0f)", hier_med, hier_med * 88.42),
          `Proximal baseline, median mg/dL (µmol/L)` = sprintf("%.2f (%.0f)", t0_med, t0_med * 88.42),
          `Median difference, mg/dL` = sprintf("%+.2f", diff_med), `AKI, hierarchy %` = f1(aki_hier),
          `AKI, proximal %` = f1(aki_t0), `Difference, pp` = sprintf("%+.1f", aki_hier - aki_t0),
          `Agreement, %` = f1(agree), `Missed by proximal, %` = f1(missed),
          `Created by proximal, %` = f1(false))], "Table4_baseline_misclassification.csv")

out(C3[, .(`Data source` = db, `Patients, n` = fn(n), `AKI-evaluable, permissive %` = f1(perm_L2),
           `AKI-evaluable, strict %` = f1(strict_L2), `AKI-evaluable, hierarchy %` = f1(hier_L2),
           `AKD-evaluable, permissive %` = f1(perm_L3), `AKD-evaluable, strict %` = f1(strict_L3),
           `AKD-evaluable, hierarchy %` = f1(hier_L3))], "Table5_ascertainment_by_rule.csv")

tim <- rbindlist(lapply(ORDER, function(k) {
  x <- TM[[k]]; if (is.null(x) || is.null(x$hour_of_day)) return(NULL)
  cc <- as.numeric(x$hour_of_day); t <- sum(cc); pk <- which.max(cc) - 1
  a <- as.numeric(x$hours_from_t0)
  data.table(`Data source` = k, `Results, n` = fn(t), `Modal hour` = sprintf("%02d:00", pk),
             `% at modal hour` = f1(100 * cc[pk + 1] / t), `00:00–05:59, %` = f1(100 * sum(cc[1:6]) / t),
             `06:00–11:59, %` = f1(100 * sum(cc[7:12]) / t), `12:00–17:59, %` = f1(100 * sum(cc[13:18]) / t),
             `18:00–23:59, %` = f1(100 * sum(cc[19:24]) / t),
             `% of window results in first 24 h after t0` = f1(100 * sum(a[25:48]) / sum(a)))
}))
out(tim, "TableS1_timing.csv")

S <- R[death_avail > 5]
out(S[, .(`Data source` = db, `Patients, n` = fn(n), `AKD-evaluable, n` = fn(L3n),
          `Died before day 7, n` = fn(died7), `Died day 7–90, n` = fn(died790),
          `Discharged, not re-tested, n` = fn(dis90), `AKD-evaluable, % of all` = f1(L3),
          `AKD-evaluable, % of measurable` = f1(akd_elig))], "TableS2_undefined_vs_unmeasured_AKD.csv")
out(S[, .(`Data source` = db, `Patients, n` = fn(n), `AKI-evaluable, %` = f1(L2),
          `Died <48 h, n` = fn(aki_died48), `Died 2–7 d, n` = fn(aki_died2_7),
          `Discharged <7 d, n` = fn(aki_disch7), `Alive, not measured, n` = fn(aki_unmeasured),
          `AKI-evaluable, % of measurable` = f1(aki_elig))], "TableS3_undefined_vs_unmeasured_AKI.csv")
cat("tables written to", normalizePath(OUT), "\n"); print(list.files(OUT))
