# Graphical abstract: 1328 x 531 px proportion (JCE), rendered at 2x.
# Every number in the text is computed from results/, so the abstract cannot
# drift from the tables when the analysis is re-run.
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(jsonlite) })
source("config.R")
dir.create("figures", showWarnings = FALSE)
R  <- fread(file.path(RESULTS_DIR, "measurement_and_provenance.csv"))
M  <- fread(file.path(RESULTS_DIR, "baseline_misclassification.csv"))
TM <- read_json(file.path(RESULTS_DIR, "timing.json"), simplifyVector = TRUE)

CC <- c("MIMIC-IV", "eICU-CRD", "AmsterdamUMCdb", "NWICU")
w <- function(d, c) sum(d[[c]] * d$n) / sum(d$n)
cc <- M[db %in% CC]; op <- M[!db %in% c(CC, "PIC")]
bars <- data.table(setting = rep(c("Critical care", "Perioperative and surgical"), each = 2),
                   baseline = rep(c("Consensus\nbaseline", "Proximal\nbaseline"), 2),
                   aki = c(w(cc, "aki_hier"), w(cc, "aki_t0"), w(op, "aki_hier"), w(op, "aki_t0")))

# ---- headline numbers
N        <- sum(R$n)
never    <- sum(R$zero * R$n) / N
has_pre  <- sum((100 - R$nopre) * R$n)
t0_share <- 100 * sum(R$t0 * R$n) / has_pre          # tier 0 among patients with any pre-index value
modal <- rbindlist(lapply(names(TM), function(k) { v <- as.numeric(TM[[k]]$hour_of_day)
  if (!length(v) || sum(v) == 0) return(NULL)
  data.table(db = k, hour = which.max(v) - 1, pct = 100 * max(v) / sum(v)) }))
top <- modal[which.max(pct)]
fmt_med <- function(x) format(x, drop0trailing = TRUE)

lines <- paste0(
  sprintf("•  Median %s to %s results per patient; %.0f%% never measured\n",
          fmt_med(min(R$med)), fmt_med(max(R$med)), never),
  sprintf("•  Drawn on early-morning blood rounds (up to %.0f%% at %02d:00)\n", top$pct, top$hour),
  sprintf("•  Consensus baseline available for %.0f%% to %.0f%% of patients\n",
          min(R$hier), max(R$hier)),
  sprintf("•  %.0f%% of patients with a pre-index value have only a <7-day one",
          t0_share))

p <- ggplot() +
  annotate("text", x = 0, y = 10.3, hjust = 0, size = 5.2, fontface = "bold", colour = "#1f3d55",
           label = "What can be measured determines what can be found") +
  annotate("text", x = 0, y = 9.2, hjust = 0, size = 3.6, colour = "grey25",
           label = sprintf("Serum creatinine in %.2f million patients, %d data sources, 6 countries",
                           N / 1e6, nrow(R))) +
  annotate("text", x = 0, y = 7.4, hjust = 0, size = 3.6, lineheight = 1.05, colour = "grey15",
           label = lines) +
  annotate("text", x = 0, y = 4.6, hjust = 0, size = 3.4, fontface = "bold", colour = "#8c2d2d",
           label = "Report measurement intensity and baseline provenance with AKI incidence") +
  coord_cartesian(xlim = c(0, 10), ylim = c(4, 11), clip = "off") + theme_void()
q <- ggplot(bars, aes(baseline, aki, fill = baseline)) + geom_col(width = .7) +
  geom_text(aes(label = sprintf("%.1f%%", aki)), vjust = -.4, size = 3) +
  facet_wrap(~setting) + scale_fill_manual(values = c("#1f3d55", "#c9a227"), guide = "none") +
  scale_y_continuous(limits = c(0, 45), expand = c(0, 0)) +
  labs(x = NULL, y = "AKI incidence, same patients (%)") + theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ragg::agg_png("figures/Graphical_abstract.png", width = 2656, height = 1062, res = 200)
print(patchwork::wrap_plots(p, q, widths = c(1.55, 1))); invisible(dev.off())
ragg::agg_tiff("figures/Graphical_abstract.tif", width = 2656, height = 1062, res = 200, compression = "lzw")
print(patchwork::wrap_plots(p, q, widths = c(1.55, 1))); invisible(dev.off())
cat("graphical abstract written\n"); cat(lines, "\n")
