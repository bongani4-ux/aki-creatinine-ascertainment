# Graphical abstract: 1328 x 531 px proportion (JCE), rendered at 2x.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
M <- fread("results/baseline_misclassification.csv")
CC <- c("MIMIC-IV", "eICU-CRD", "AmsterdamUMCdb", "NWICU")
w <- function(d, c) sum(d[[c]] * d$n) / sum(d$n)
cc <- M[db %in% CC]; op <- M[!db %in% c(CC, "PIC")]
bars <- data.table(setting = rep(c("Critical care", "Perioperative and surgical"), each = 2),
                   baseline = rep(c("Consensus\nbaseline", "Proximal\nbaseline"), 2),
                   aki = c(w(cc, "aki_hier"), w(cc, "aki_t0"), w(op, "aki_hier"), w(op, "aki_t0")))
p <- ggplot() +
  annotate("text", x = 0, y = 10.3, hjust = 0, size = 5.2, fontface = "bold", colour = "#1f3d55",
           label = "What can be measured determines what can be found") +
  annotate("text", x = 0, y = 9.2, hjust = 0, size = 3.6, colour = "grey25",
           label = "Serum creatinine in 1.26 million patients, 14 data sources, 6 countries") +
  annotate("text", x = 0, y = 7.4, hjust = 0, size = 3.6, lineheight = 1.05, colour = "grey15",
           label = paste0("\u2022  Measured 0 to 12 times per patient; 26% never measured\n",
                          "\u2022  Drawn on early-morning blood rounds (up to 74% at 06:00)\n",
                          "\u2022  Consensus baseline available for 0% to 74% of patients\n",
                          "\u2022  Nearly half of baselines come from an unspecified window")) +
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
cat("graphical abstract written\n")
