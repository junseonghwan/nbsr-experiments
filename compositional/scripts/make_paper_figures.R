#!/usr/bin/env Rscript
# Rebuild the paper's benchmark figures from the released per-feature outputs (the Zenodo deposit).
#
#   Rscript compositional/scripts/make_paper_figures.R [<data_dir>] [<out_dir>]
#   defaults: compositional/data/outputs (the unpacked deposit) and compositional/figures
#
# Inputs (one row per method x case x n x replicate x feature, unfiltered; see the deposit's README):
#   mirna_{nbsr,compositional,countbased}_features.csv.gz   miRNA benchmark: NBSR (HMC, Laplace); LinDA, ALDEx2, ANCOM-BC(2); DESeq2, edgeR
#   microbiome_features.csv.gz   microbiome benchmark (LNDM), NBSR (Laplace, diffuse intercept prior) and the competitors
#   conc_<A>_features.csv.gz     the concentration sweep (plain Dirichlet-multinomial at A = 38 ... 4000, and 10x depth)
# Outputs (images/):
#   fig_ranking_<bench>_pr.pdf, fig_ranking_<bench>_ap.pdf   precision-recall and mean average precision
#   fig_qvalue_<bench>.png, fig_lfsr_<bench>.png              adjusted p / lfsr against the true fold change
#   fig2a_coverage_vs_nominal.pdf, fig2b_zsd_by_abundance.pdf, fig3_interval_geometry.pdf   calibration (miRNA, HMC)
#   fig_concentration.pdf                                       samples needed vs concentration, and the sweep
# Not rebuilt here: Figure 2 (fig_rare_feature) needs the posterior draws of two fits, which are not in the deposit,
# and the compositional-shift toy figure is drawn, not computed.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "compositional/data/outputs"
out <- if (length(args) >= 2) args[2] else "compositional/figures"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
blue <- "#2a78d6"; orange <- "#eb6834"; grey_m <- "#8a8985"; grey_l <- "#c3c2bb"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e0"
labels <- c("nbsr-hmc" = "NBSR (HMC)", "nbsr-laplace" = "NBSR (Laplace)", "nbsr-laplace-wideint" = "NBSR (Laplace)",
            "linda" = "LinDA", "aldex2" = "ALDEx2", "ancombc" = "ANCOM-BC", "ancombc2" = "ANCOM-BC2",
            "deseq2" = "DESeq2", "edger-qlf" = "edgeR QLF", "edger-lrt" = "edgeR LRT")
cols <- c("NBSR (HMC)" = blue, "NBSR (Laplace)" = orange, "DESeq2" = grey_m, "edgeR QLF" = grey_m, "edgeR LRT" = grey_m,
          "LinDA" = grey_l, "ANCOM-BC" = grey_l, "ANCOM-BC2" = grey_l, "ALDEx2" = grey_l)
theme_paper <- function(base_size = 8) theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
        axis.title = element_text(colour = ink2), axis.text = element_text(colour = ink), strip.text = element_text(colour = ink, face = "bold"),
        legend.position = "none", plot.title = element_text(face = "bold", size = rel(1), colour = ink), plot.title.position = "plot",
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 6, 4, 4))
pct <- function(x) paste0(round(100 * x), "%")
KEY <- c("simul_case", "sample_size", "replicate_id", "feature"); REP <- c("method", "simul_case", "sample_size", "replicate_id")
cfgs <- list(
  mirna = list(lfsr_method = "nbsr-hmc", common = FALSE, cases = c("LOW", "MED", "HIGH"),
               methods = c("nbsr-hmc", "nbsr-laplace", "linda", "aldex2", "ancombc", "ancombc2", "deseq2", "edger-qlf", "edger-lrt")),
  microbiome = list(lfsr_method = "nbsr-laplace-wideint", common = TRUE, cases = c("ALL25", "ALL50", "BAL50"),
               methods = c("nbsr-laplace-wideint", "linda", "aldex2", "ancombc2")))

read_gz <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)))   # no R.utils dependency
read_features <- function(bench) {   # the miRNA benchmark is released in three files (GitHub's file-size limit)
  fs <- Sys.glob(file.path(data_dir, paste0(bench, c("", "_nbsr", "_compositional", "_countbased"), "_features.csv.gz")))
  rbindlist(lapply(fs, read_gz), fill = TRUE) }
load_bench <- function(bench) {
  cfg <- cfgs[[bench]]
  D <- read_features(bench)[method %in% cfg$methods & simul_case %in% cfg$cases]
  for (col in c("tested", "called", "passed_sens")) if (!col %in% names(D)) D[, (col) := NA]
  for (col in c("q_value", "estimate", "se", "p_value", "true_effect_abs")) if (!col %in% names(D)) D[, (col) := NA_real_]
  nb <- startsWith(D$method, "nbsr")
  D[nb, `:=`(tested = TRUE, estimate = post_mean, se = post_sd, p_value = 2 * pnorm(-abs(post_mean / post_sd)), lfsr = pmin(p_pos, 1 - p_pos))]
  D[!nb, lfsr := NA_real_]
  D[is.na(tested), tested := FALSE]; D[is.na(true_effect_abs), true_effect_abs := true_effect]
  D[nb & tested == TRUE, q_value := p.adjust(p_value, "BH"), by = REP]          # NBSR's adjusted p: posterior-z Wald + BH
  D[, call := fcase(startsWith(method, "nbsr"), tested & !is.na(q_value) & q_value <= 0.05,
                    method == "ancombc2", tested & called %in% TRUE & passed_sens %in% TRUE, default = tested & called %in% TRUE)]
  D[is.na(call), call := FALSE]
  if (cfg$common) { ok <- D[, .(ok = all(tested) && .N == length(cfg$methods)), by = KEY][ok == TRUE, ..KEY]; D <- merge(D, ok, by = KEY) }
  D[, stat := fifelse(tested & is.finite(estimate / se), abs(estimate / se), -Inf)]
  D[, method_lab := factor(labels[method], labels[cfg$methods])]
  cat(bench, ":", nrow(D), "rows,", paste(cfg$methods, collapse = " "), "\n"); D
}

ranking_figures <- function(bench) {
  cfg <- cfgs[[bench]]; D <- load_bench(bench); ns <- sort(unique(D$sample_size))
  pr_curve <- function(stat, de) { o <- order(-stat); de <- de[o]; tp <- cumsum(de); k <- seq_along(de)
    i <- unique(c(seq(1, length(de), length.out = min(length(de), 600)), length(de))); data.table(recall = tp[i] / sum(de), precision = tp[i] / k[i]) }
  PR <- D[, pr_curve(stat, is_de), by = .(method_lab, sample_size)]
  ap_one <- function(stat, de) { o <- order(-stat); de <- de[o]; tp <- cumsum(de); sum((tp / seq_along(de))[de]) / sum(de) }
  AP <- D[, .(ap = ap_one(stat, is_de)), by = c("method_lab", REP)][, .(ap = mean(ap)), by = .(method_lab, sample_size)]
  is_nbsr <- function(m) grepl("^NBSR", m)
  PR <- rbind(PR[!is_nbsr(method_lab)], PR[is_nbsr(method_lab)])
  p_a <- ggplot(PR, aes(recall, precision, colour = method_lab, group = method_lab)) + geom_path(aes(linewidth = is_nbsr(method_lab))) +
    facet_wrap(~sample_size, nrow = 1, labeller = labeller(sample_size = function(x) paste0("n = ", x))) +
    scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.9)) + scale_colour_manual(values = cols) +
    scale_x_continuous(labels = pct, limits = c(0, 1)) + scale_y_continuous(labels = pct, limits = c(0, 1)) +
    labs(title = "Precision against recall, pooled over replicates", x = "recall (perturbed features found)", y = "precision") + theme_paper()
  lab <- AP[sample_size == max(ns)][order(ap)]; lab[, y := ap]; gap <- 0.045
  for (i in seq_len(nrow(lab))[-1]) if (lab$y[i] - lab$y[i - 1] < gap) lab$y[i] <- lab$y[i - 1] + gap
  lab[, y := y - max(0, max(y) - 1)]; lab[, colour := ifelse(is_nbsr(method_lab), cols[as.character(method_lab)], ink2)]
  xr <- range(ns); xpad <- diff(xr) * 0.55
  p_b <- ggplot(AP, aes(sample_size, ap, colour = method_lab, group = method_lab)) + geom_line(aes(linewidth = is_nbsr(method_lab))) + geom_point(aes(size = is_nbsr(method_lab))) +
    geom_segment(data = lab, aes(x = max(ns) * 1.01, xend = max(ns) * 1.06, y = ap, yend = y), colour = grid_col, linewidth = 0.3) +
    geom_text(data = lab, aes(x = max(ns) * 1.07, y = y, label = method_lab), colour = lab$colour, hjust = 0, size = 2.3) +
    scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.9)) + scale_size_manual(values = c(`FALSE` = 0.9, `TRUE` = 1.6)) + scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = ns, limits = c(xr[1] * 0.9, xr[2] + xpad)) + scale_y_continuous(labels = pct, limits = c(0, 1)) + coord_cartesian(clip = "off") +
    labs(title = "Mean average precision", x = "samples per condition", y = "average precision") + theme_paper()
  ggsave(file.path(out, paste0("fig_ranking_", bench, "_pr.pdf")), p_a, width = 1.9 * length(ns) + 0.6, height = 2.4, device = cairo_pdf)
  ggsave(file.path(out, paste0("fig_ranking_", bench, "_ap.pdf")), p_b, width = 3.0, height = 2.4, device = cairo_pdf)
  cat("mean AP:\n"); print(dcast(AP[, .(method_lab, sample_size, v = sprintf("%.3f", ap))], method_lab ~ sample_size, value.var = "v"))

  outcome_cols <- c(TP = "#5aa469", FN = "#eb9a5c", FP = "#c0392b", TN = "#5b8fd1")
  Q <- D[tested == TRUE & !is.na(q_value)]
  Q[, outcome := factor(fcase(call & is_de, "TP", call & !is_de, "FP", !call & is_de, "FN", default = "TN"), c("TP", "FN", "FP", "TN"))]
  Q[, `:=`(q_plot = pmax(q_value, 1e-4), x = true_effect_abs / log(2))]
  Q <- Q[order(match(outcome, c("TN", "FN", "TP", "FP")))]
  st <- Q[, .(label = sprintf("FDR %.2f\nTPR %.2f", sum(call & !is_de) / max(sum(call), 1), sum(call & is_de) / sum(is_de))), by = .(method_lab, sample_size)]
  p_q <- ggplot(Q, aes(x, q_plot, colour = outcome)) + geom_hline(yintercept = 0.05, colour = "#c0392b", linetype = "22", linewidth = 0.35) +
    geom_point(size = 0.35, alpha = 0.6, position = position_jitter(width = 0.03, height = 0, seed = 1)) +
    geom_text(data = st, aes(x = Inf, y = 1, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.1, size = 2.1, colour = ink, lineheight = 0.9) +
    scale_colour_manual(values = outcome_cols, labels = c(TP = "true positive", FN = "false negative", FP = "false positive", TN = "true negative")) +
    scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2, 0.05, 1), labels = c(expression("" <= 10^-4), expression(10^-3), expression(10^-2), "0.05", "1")) +
    facet_grid(method_lab ~ sample_size, labeller = labeller(sample_size = function(x) paste0("n = ", x))) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(x = expression("true " * log[2] * " fold change (absolute scale)"), y = "adjusted p-value") +
    theme_paper() + theme(legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.key.size = unit(8, "pt"), strip.text.y = element_text(angle = 0))
  ggsave(file.path(out, paste0("fig_qvalue_", bench, ".png")), p_q, width = 1.6 * length(ns) + 1.3, height = 0.95 * length(cfg$methods) + 0.9, dpi = 300, device = ragg::agg_png, bg = "white")

  L <- D[method == cfg$lfsr_method]; floor_lfsr <- 1 / 16000; thr <- 1e-3
  L[, call_l := pmax(p_pos, 1 - p_pos) >= 1 - thr]
  L[, outcome := factor(fcase(call_l & is_de, "TP", call_l & !is_de, "FP", !call_l & is_de, "FN", default = "TN"), c("TP", "FN", "FP", "TN"))]
  L[, `:=`(lfsr_plot = pmax(lfsr, floor_lfsr), x = true_effect_abs / log(2), setting = factor(simul_case, cfg$cases))]
  L <- L[order(match(outcome, c("TN", "FN", "TP", "FP")))]
  stl <- L[, .(label = sprintf("FDR %.1f%%\nFNR %.1f%%", 100 * sum(call_l & !is_de) / max(sum(call_l), 1), 100 * sum(!call_l & is_de) / sum(is_de))), by = .(setting, sample_size)]
  p_l <- ggplot(L, aes(x, lfsr_plot, colour = outcome)) + geom_hline(yintercept = thr, colour = "#c0392b", linetype = "22", linewidth = 0.35) +
    geom_point(size = 0.4, alpha = 0.6, position = position_jitter(width = 0.03, height = 0, seed = 1)) +
    geom_text(data = stl, aes(x = Inf, y = 0.45, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1, size = 2.1, colour = ink, lineheight = 0.9) +
    scale_colour_manual(values = outcome_cols, labels = c(TP = "true positive", FN = "false negative", FP = "false positive", TN = "true negative")) +
    scale_y_log10(breaks = c(floor_lfsr, 1e-3, 1e-2, 1e-1, 0.5), labels = c("0", expression(10^-3), expression(10^-2), "0.1", "0.5")) +
    facet_grid(sample_size ~ setting, scales = "free_x", labeller = labeller(sample_size = function(x) paste0("n = ", x))) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(x = expression("true " * log[2] * " fold change (absolute scale)"), y = "local false sign rate") +
    theme_paper() + theme(legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.key.size = unit(8, "pt"), strip.text.y = element_text(angle = 0))
  ggsave(file.path(out, paste0("fig_lfsr_", bench, ".png")), p_l, width = 6.5, height = 1.5 * length(ns) + 0.8, dpi = 300, device = ragg::agg_png, bg = "white")
}
ranking_figures("mirna"); ranking_figures("microbiome")

## ---- calibration figures (miRNA, NBSR by HMC): report_findings.R ---------------------------------------------------
d <- read_features("mirna")[method == "nbsr-hmc"]
d[, label := factor("nbsr-hmc")]
d[, err := post_mean - true_effect]; d[, z := err / post_sd]; d[, w95 := q975 - q025]
d[, c50 := q25 <= true_effect & true_effect <= q75]; d[, c90 := q05 <= true_effect & true_effect <= q95]
d[, c95 := q025 <= true_effect & true_effect <= q975]; d[, c99 := q005 <= true_effect & true_effect <= q995]
d[, call95 := q025 > 0 | q975 < 0]; d[, edge := pmin(abs(q025), abs(q975))]
d[, stratum := fifelse(is_de, "DE", "non-DE")]
by_unit <- c("label", "sample_size")
d[is_de == FALSE, abund := cut(alpha_null, quantile(alpha_null, 0:5 / 5), include.lowest = TRUE,
                              labels = c("Q1 lowest", "Q2", "Q3", "Q4", "Q5 highest")), by = by_unit]
SLOTS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
label_cols <- setNames(SLOTS[1], "nbsr-hmc")
pair_cols <- c("non-DE" = SLOTS[1], "DE" = SLOTS[2], "false positive" = SLOTS[1], "true positive" = SLOTS[2])
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
  axis.line = element_blank(), strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom", plot.title.position = "plot"))
size_lab <- function(n) paste0("n = ", n, " per group")
cov_long <- melt(d[, .(cov50 = mean(c50), cov90 = mean(c90), cov95 = mean(c95), cov99 = mean(c99)), by = c(by_unit, "stratum")],
                 id.vars = c(by_unit, "stratum"), variable.name = "level", value.name = "coverage")
cov_long[, nominal := as.numeric(sub("cov", "", level)) / 100]
p2a <- ggplot(cov_long, aes(nominal, coverage, colour = stratum)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9a9994") + geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_grid(label ~ sample_size, labeller = labeller(sample_size = size_lab)) +
  scale_colour_manual(values = pair_cols[c("non-DE", "DE")], name = NULL) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(title = "Empirical coverage against nominal level", subtitle = "above the dashed line = conservative", x = "nominal", y = "empirical")
zq <- d[is_de == FALSE, .(z_sd = sd(z)), by = c(by_unit, "abund")]
p2b <- ggplot(zq, aes(abund, z_sd, colour = label, group = label)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#9a9994") + geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_wrap(~ sample_size, labeller = as_labeller(size_lab)) + scale_colour_manual(values = label_cols, name = NULL) +
  labs(title = "Posterior width relative to actual error, non-DE features",
       subtitle = "sd of (post_mean - truth) / post_sd; 1 = calibrated, below 1 = conservative", x = "abundance quintile (alpha_null)", y = "z-score sd")
ggsave(file.path(out, "fig2a_coverage_vs_nominal.pdf"), p2a, width = 7.5, height = 1.6 + 2.2)
ggsave(file.path(out, "fig2b_zsd_by_abundance.pdf"), p2b, width = 7, height = 4)
d[, outcome := fifelse(is_de, "true positive", "false positive")]
p3 <- ggplot(d[call95 == TRUE], aes(edge / w95, w95, colour = outcome)) + geom_point(size = 1.4, alpha = 0.55) +
  facet_grid(label ~ sample_size, labeller = labeller(sample_size = size_lab)) +
  scale_colour_manual(values = pair_cols[c("false positive", "true positive")], name = NULL) + scale_x_log10() + scale_y_log10() +
  labs(title = "Discoveries at the 95% interval: how far past zero, and how wide",
       x = "distance of nearest interval edge from 0, as a fraction of interval width (log)", y = "95% interval width (log)")
ggsave(file.path(out, "fig3_interval_geometry.pdf"), p3, width = 7, height = 2.5 + 2)

## ---- concentration figure: concentration_figure.R --------------------------------------------------------------------
z <- qnorm(0.975) + qnorm(0.8)
n_req <- function(pi, s, A) { I <- pi * s * A / (s + A); 2 * z^2 / (I * log(2)^2) }
grid <- 10^seq(-5, -0.7, length.out = 200)
curves <- rbindlist(list(
  data.table(pi = grid, n = n_req(grid, 9000, 38),   setting = "microbiome template: A = 38, 9,000 reads"),
  data.table(pi = grid, n = n_req(grid, 9000, 100),  setting = "benchmark: A = 100, 9,000 reads"),
  data.table(pi = grid, n = n_req(grid, 90000, 100), setting = "benchmark at 10x depth: A = 100, 90,000 reads"),
  data.table(pi = grid, n = n_req(grid, 1e7, 4350),  setting = "miRNA benchmark: A = 4,350, 10^7 reads")))
curves[, setting := factor(setting, unique(setting))]
pa <- ggplot(curves, aes(pi, n, colour = setting, linetype = setting)) +
  annotate("rect", xmin = 1e-5, xmax = 1e-3, ymin = 1, ymax = 3e5, alpha = 0.06) +
  annotate("text", x = c(10^-4, 10^-2.5, 10^-1.35), y = 1.6e5, label = c("LOW", "MED", "HIGH"), size = 3, colour = "grey30") +
  geom_vline(xintercept = c(1e-3, 1e-2), colour = "grey80", linewidth = 0.4) + geom_line(linewidth = 0.9) +
  scale_x_log10(breaks = 10^(-5:-1), labels = c("0.001%", "0.01%", "0.1%", "1%", "10%")) +
  scale_y_log10(breaks = 10^(0:5), labels = c("1", "10", "100", "1,000", "10,000", "100,000")) + coord_cartesian(ylim = c(1, 3e5)) +
  scale_colour_manual(values = c("#c4502a", "#2f6f7e", "#8fc3cb", "#5b5b5b")) + scale_linetype_manual(values = c("solid", "solid", "22", "solid")) +
  labs(x = "mean abundance of the feature", y = "samples per group for a two-fold change", title = "(a) Samples needed for a two-fold change") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "vertical",
                                        panel.grid.minor = element_blank(), plot.title = element_text(size = 10, face = "bold"))
gens <- c("conc_A38" = "A = 38", "conc_A100" = "A = 100", "conc_A400" = "A = 400", "conc_A1600" = "A = 1,600", "conc_A4000" = "A = 4,000", "conc_A100_d10" = "A = 100, 10x depth")
pw <- rbindlist(lapply(names(gens), function(g) {
  f <- file.path(data_dir, paste0(g, "_features.csv.gz")); if (!file.exists(f)) return(NULL)
  s <- read_gz(f)[method == "nbsr-laplace-wideint" & simul_case == "ALL25" & is_de == TRUE]
  s[, .(TPR = mean(pmax(p_pos, 1 - p_pos) >= 0.999)), by = .(sample_size, stratum)][, setting := gens[[g]]] }))
pw[, `:=`(setting = factor(setting, gens), stratum = factor(stratum, c("HIGH", "MED", "LOW")))]
pb <- ggplot(pw, aes(sample_size, TPR, colour = setting, linetype = setting)) + geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  facet_wrap(~stratum, nrow = 1) + scale_x_log10(breaks = c(40, 80, 160, 320)) + scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c("#c4502a", "#2f6f7e", "#6a8f6f", "#8a6fa8", "#5b5b5b", "#8fc3cb")) + scale_linetype_manual(values = c(rep("solid", 5), "22")) +
  labs(x = "samples per group", y = "power (NBSR, lfsr <= 0.001)", title = "(b) Power rises with the concentration far more than with depth") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank(),
                                        plot.title = element_text(size = 10, face = "bold")) + guides(colour = guide_legend(nrow = 1))
pdf(file.path(out, "fig_concentration.pdf"), width = 10, height = 4.8)
grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(0.42, 0.58))))
print(pa, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)); print(pb, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())
cat("wrote the figures to", out, "\n")
