# =============================================================================
# figures.R  -  build all manuscript figures from results/.
# Usage:  Rscript figures.R
# Writes figures/Figure1-5 and SuppFigure1, each as vector PDF (preferred by
# the journal for charts) and as 600 dpi LZW-compressed TIFF.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(ggplot2); library(patchwork); library(scales)
})
source("config.R")
OUT <- "figures"; dir.create(OUT, showWarnings = FALSE)

COL <- c("Critical care" = "#31506e", "Critical care (paed)" = "#7fa3c4",
         "Perioperative" = "#b8562f", "Surgical (inst)" = "#6b3a63")
theme_pub <- function(base = 8) theme_classic(base_size = base, base_family = "sans") +
  theme(panel.grid.major.x = element_line(colour = "grey88", linewidth = .3),
        axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .3),
        plot.title = element_text(size = base + .5, face = "plain", hjust = 0),
        plot.title.position = "plot", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = base - 1.5), legend.title = element_blank(),
        legend.background = element_blank())

save_fig <- function(p, name, w = 7.5, h = 4.6) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ragg::agg_tiff(file.path(OUT, paste0(name, ".tif")), width = w, height = h, units = "in",
                 res = 600, compression = "lzw"); print(p); invisible(dev.off())
}

R  <- fread(file.path(RESULTS_DIR, "measurement_and_provenance.csv"))
C3 <- fread(file.path(RESULTS_DIR, "ascertainment_by_baseline_rule.csv"))
M  <- fread(file.path(RESULTS_DIR, "baseline_misclassification.csv"))
TM <- read_json(file.path(RESULTS_DIR, "timing.json"), simplifyVector = TRUE)
SET <- setNames(R$setting, R$db)

# ---------------------------------------------------------------- Figure 1
d <- copy(R)[order(med, db)]; d[, db := factor(db, levels = db)]
pA <- ggplot(d, aes(y = db, colour = setting)) +
  geom_segment(aes(x = q1, xend = q3, yend = db), linewidth = 2.6, alpha = .35) +
  geom_segment(aes(x = q3, xend = p95, yend = db), linewidth = .5, alpha = .7) +
  geom_point(aes(x = med), size = 1.9) +
  geom_text(aes(x = p95 + 2, label = ifelse(abs(med - round(med)) < 1e-9, sprintf("%d", as.integer(round(med))), sprintf("%.1f", med))), hjust = 0, size = 2.2, colour = "grey30") +
  scale_colour_manual(values = COL, guide = "none") + scale_x_continuous(limits = c(0, 78), expand = c(0, 0)) +
  labs(x = "Creatinine results per patient\n(median, IQR, 95th centile)", y = NULL,
       title = "(A) Measurement intensity") + theme_pub()
d2 <- copy(R)[order(zero, db)]; d2[, db := factor(db, levels = db)]
pB <- ggplot(d2, aes(x = zero, y = db, fill = setting)) + geom_col(width = .66) +
  geom_text(aes(label = sprintf("%.1f", zero)), hjust = -.15, size = 2.2) +
  scale_fill_manual(values = COL) + scale_x_continuous(limits = c(0, 72), expand = c(0, 0)) +
  labs(x = "Patients with no creatinine measured\nfrom 365 d before to 90 d after t\u2080 (%)", y = NULL,
       title = "(B) Never measured") + theme_pub() + theme(legend.position = c(.72, .18))
save_fig(pA + pB, "Figure1")

# ---------------------------------------------------------------- Figure 2
hc <- c("HiRID" = "#1f3d55", "NWICU" = "#4a7fa5", "MOVER" = "#b8562f", "MIMIC-IV" = "#31506e", "PIC" = "#7fa3c4")
hod <- rbindlist(lapply(names(hc), function(k) if (!is.null(TM[[k]]$hour_of_day)) {
  v <- as.numeric(TM[[k]]$hour_of_day); data.table(db = k, hour = 0:23, pct = 100 * v / sum(v)) }))
hod[, db := factor(db, levels = names(hc))]
p2a <- ggplot(hod, aes(hour, pct, colour = db)) +
  annotate("rect", xmin = 4, xmax = 7, ymin = -Inf, ymax = Inf, fill = "#c9a227", alpha = .15) +
  annotate("text", x = 5.5, y = max(hod$pct) * .82, label = "ward round\nblood draw", size = 2.1, colour = "#8a6d00") +
  geom_line(linewidth = .55) + geom_point(size = .9) + scale_colour_manual(values = hc) +
  scale_x_continuous(breaks = seq(0, 21, 3), labels = sprintf("%02d", seq(0, 21, 3))) +
  labs(x = "Hour of day (clock time)", y = "Percentage of creatinine results",
       title = "(A) When in the day is creatinine drawn?") + theme_pub() +
  theme(legend.position = c(.82, .75), panel.grid.major.y = element_line(colour = "grey88", linewidth = .3))
fc <- c("MIMIC-IV" = "#31506e", "eICU-CRD" = "#5a86ad", "HiRID" = "#1f3d55", "SICdb" = "#9ab4c9",
        "INSPIRE" = "#c9773f", "MOVER" = "#b8562f")
hft <- rbindlist(lapply(names(fc), function(k) { a <- as.numeric(TM[[k]]$hours_from_t0)
  data.table(db = k, h = seq(-24, 167) + .5, pct = 100 * a / max(sum(a), 1)) }))
hft[, db := factor(db, levels = names(fc))]
p2b <- ggplot(hft, aes(h, pct, colour = db)) + geom_vline(xintercept = 0, linewidth = .35, colour = "grey20") +
  geom_line(linewidth = .45) + scale_colour_manual(values = fc) +
  scale_x_continuous(breaks = seq(-24, 168, 24), limits = c(-24, 168)) +
  labs(x = expression("Hours from t"[0]), y = "Percentage of results in window, per hour",
       title = "(B) Hourly timing around the index episode") + theme_pub() +
  theme(legend.position = c(.82, .75), panel.grid.major.y = element_line(colour = "grey88", linewidth = .3))
save_fig(p2a + p2b, "Figure2", h = 3.9)

# ---------------------------------------------------------------- Figure 3
d <- copy(R)[order(-hier, db)]; lev <- d$db
long <- melt(d[, .(db, T1, T2, T3, T4, t0, nopre)], id.vars = "db")
tl <- c(T1 = "Tier 1 (7\u201330 d)", T2 = "Tier 2 (31\u201390 d)", T3 = "Tier 3 (91\u2013180 d)",
        T4 = "Tier 4 (181\u2013365 d)", t0 = "Tier 0 (<7 d, unspecified by consensus)", nopre = "No pre-index value")
tc <- c("#1f3d55", "#3d6f96", "#7ba0bf", "#b6cbdb", "#c9a227", "#dcdcdc"); names(tc) <- tl
long[, `:=`(db = factor(db, levels = lev), variable = factor(tl[as.character(variable)], levels = rev(tl)))]
p3 <- ggplot(long, aes(db, value, fill = variable)) + geom_col(width = .72) +
  geom_text(data = d, aes(x = factor(db, levels = lev), y = 104, label = sprintf("%.0f", hier)),
            inherit.aes = FALSE, size = 2.3, fontface = "bold", colour = "#1f3d55") +
  annotate("text", x = .4, y = 110, label = "consensus-compliant %:", hjust = 0, size = 2.2,
           fontface = "bold", colour = "#1f3d55") +
  scale_fill_manual(values = tc, breaks = tl) + scale_y_continuous(limits = c(0, 112), expand = c(0, 0)) +
  labs(x = NULL, y = "Percentage of patients") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = .3), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
save_fig(p3, "Figure3", w = 7.2, h = 4.6)

# ---------------------------------------------------------------- Figure 4
d <- copy(M)[order(aki_hier, db)]; d[, `:=`(db = factor(db, levels = db), setting = SET[as.character(db)])]
p4a <- ggplot(d, aes(y = db)) +
  geom_segment(aes(x = aki_t0, xend = aki_hier, yend = db), colour = "grey60", linewidth = .4) +
  geom_point(aes(x = aki_t0, colour = setting), shape = 21, fill = "white", size = 1.9, stroke = .7) +
  geom_point(aes(x = aki_hier, colour = setting), size = 1.9) +
  geom_text(aes(x = pmax(aki_hier, aki_t0) + 1.4, label = sprintf("%+.1f", aki_hier - aki_t0)),
            hjust = 0, size = 2.2, colour = "grey30") +
  scale_colour_manual(values = COL, guide = "none") + scale_x_continuous(limits = c(0, 52), expand = c(0, 0)) +
  labs(x = "AKI incidence in the same patients (%)", y = NULL, title = "(A) Effect of baseline choice",
       caption = "Filled: hierarchy baseline   Open: proximal (tier 0) baseline") +
  theme_pub() + theme(plot.caption = element_text(size = 6, hjust = 1))
mm <- melt(d[, .(db, missed, false)], id.vars = "db")
mm[, variable := factor(variable, levels = c("missed", "false"),
       labels = c("AKI missed by proximal baseline", "AKI created by proximal baseline"))]
p4b <- ggplot(mm, aes(value, db, fill = variable)) +
  geom_col(position = position_dodge2(reverse = TRUE, padding = .05), width = .78) +
  geom_text(data = d, aes(x = missed + .4, y = db, label = sprintf("%.1f", missed)), inherit.aes = FALSE,
            hjust = 0, size = 2.2, nudge_y = .19) +
  scale_fill_manual(values = c("#8c2d2d", "#c9a227")) + scale_x_continuous(limits = c(0, 22), expand = c(0, 0)) +
  labs(x = "Patients misclassified (%)", y = NULL, title = "(B) Direction of misclassification") +
  theme_pub() + theme(legend.position = c(.66, .1))
save_fig(p4a + p4b, "Figure4", h = 4.4)

# ---------------------------------------------------------------- Figure 5
d <- copy(C3)[order(strict_L2, db)]; lev <- d$db
rc <- c(perm = "#c9a227", strict = "#5a86ad", hier = "#1f3d55")
rl <- c(perm = "Permissive: any value up to 24 h after t\u2080", strict = "Strict: any value before t\u2080",
        hier = "Consensus hierarchy: 7\u2013365 d before t\u2080")
mk <- function(suf, ttl, xmax, legend) {
  x <- melt(d[, c("db", paste0(c("perm", "strict", "hier"), "_", suf)), with = FALSE], id.vars = "db")
  # plotting order hier, strict, perm puts the permissive bar on top within each source
  x[, `:=`(rule = factor(sub(paste0("_", suf), "", variable), levels = c("hier", "strict", "perm")),
           db = factor(db, levels = lev))]
  dodge <- position_dodge(width = .8)
  ggplot(x, aes(value, db, fill = rule, group = rule)) +
    geom_col(position = dodge, width = .78) +
    geom_text(data = x[rule != "strict"], aes(label = sprintf("%.0f", value), colour = rule),
              position = dodge, hjust = -.3, size = 1.9, show.legend = FALSE) +
    scale_fill_manual(values = rc, labels = rl, breaks = c("perm", "strict", "hier")) +
    scale_colour_manual(values = c(perm = "#8a6d00", strict = "#5a86ad", hier = "#1f3d55")) +
    scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
    labs(x = "Percentage of patients", y = NULL, title = ttl) + theme_pub() +
    theme(legend.position = if (legend) "bottom" else "none")
}
p5 <- mk("L2", "(A) AKI-evaluable", 112, TRUE) +
  (mk("L3", "(B) AKD-evaluable", 60, FALSE) + theme(axis.text.y = element_blank())) +
  plot_layout(widths = c(1.7, 1), guides = "collect") & theme(legend.position = "bottom")
save_fig(p5, "Figure5", h = 5.2)

# ---------------------------------------------------------------- Supplementary Figure 1
sub <- copy(R)[as.logical(death_recorded) == TRUE & L1n > 0][order(L3, db)]
cats <- c(L3n = "AKD-evaluable", nopre_n = "No pre-index value",
          died7 = "Died before day 7 (undefined)", died790 = "Died day 7–90",
          dis90 = "Discharged, not re-tested", inh = "Other / not classifiable")
fills <- setNames(c("#1f3d55", "#dcdcdc", "#6e6e6e", "#a5a5a5", "#8c2d2d", "#c96a6a"), cats)
seg <- melt(sub[, .(db, n, L3n, nopre_n, died7, died790, dis90, inh)], id.vars = c("db", "n"))
# reversed levels so the first category (AKD-evaluable) stacks at the left
seg[, `:=`(pct = 100 * value / n, db = factor(db, levels = sub$db),
           variable = factor(cats[as.character(variable)], levels = rev(cats)))]
ps <- ggplot(seg, aes(pct, db, fill = variable)) + geom_col(width = .62) +
  geom_text(data = sub, aes(x = 101, y = factor(db, levels = sub$db), label = sprintf("%.1f%%", akd_elig)),
            inherit.aes = FALSE, hjust = 0, size = 2.3, fontface = "bold", colour = "#1f3d55") +
  scale_fill_manual(values = fills, breaks = cats) +
  scale_x_continuous(limits = c(0, 112), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "Percentage of patients", y = NULL,
       caption = "Bold: AKD-evaluable among patients with a pre-index value who were alive at day 90") +
  theme_pub() + theme(legend.position = "bottom", plot.caption = element_text(size = 6)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
save_fig(ps, "SuppFigure1", w = 7, h = 4)
cat("figures written to", normalizePath(OUT), "\n")
