#!/usr/bin/env Rscript
# Rebuild the paper's tables from the released per-feature outputs (the Zenodo deposit).
#
#   Rscript recomb2026/scripts/paper_tables.R [<data_dir>] [<out_dir>]
#   defaults: recomb2026/data/outputs (the unpacked deposit) and recomb2026/tables
#
# Every table that the paper computes from per-feature summaries is written as a CSV and collected in tables.md:
#   table1_mirna.csv         Table 1: TPR / FDR at each method's own rule and 95% coverage of the perturbed
#                            features, miRNA benchmark; plus median average precision
#   table2_microbiome.csv    Table 2 and the appendix grid: power / FDR on the genera every method tested, each
#                            method at its stated rule, and NBSR's lfsr and delta-method rows
#   matched_fdr.csv          power at realised FDR 0.05 and 0.10 on the common set (appendix)
#   dispersion_variants.csv  the per-sample dispersion models on the fixed- and drawn-concentration generators
#   table3_selex.csv         Table 3: calls on the 27 characterised variants and by direction on the rest
#   permutation.csv          the label-permutation null: false positives per replicate, share of replicates with
#                            a call, raw p below 0.05
#   real_labels.csv          the case-control studies with their real labels: calls by direction, the
#                            established colorectal-cancer genera
#   calibration.csv          appendix: coverage, z-sd and RMSE of NBSR (HMC) by n and feature type
#   decision_rules.csv       appendix: the interval, lfsr and Bayesian-FDR rules on identical posteriors
#   laplace_vs_hmc.csv       appendix: sampler against Laplace approximation (and its empirical-Bayes prior)
#   wald_rules.csv           appendix: delta-method, posterior-z and lfsr rules with BH
#   ablations.csv            appendix: every ablation variant by n (coverage, bias, calls, ranking)
# Not reproducible from summaries (they need the posterior draws): the variance decomposition of the shift
# (tab:gap), the reference-set sensitivity (tab:cutoff, tab:abund), the dispersion-regularisation factorial on
# the fitted composition (tab:regularise) and Figure 2.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "recomb2026/data/outputs"
out <- if (length(args) >= 2) args[2] else "recomb2026/tables"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ALPHA <- 0.05; LFSR_THR <- 1e-3; LFSR_TOL <- 1e-9
read_gz <- function(f) fread(cmd = paste("gzip -dc", shQuote(f)))
read_set <- function(stem) {   # a dataset released as one or several files
  fs <- Sys.glob(file.path(data_dir, paste0(stem, c("", "_*"), "_features.csv.gz")))
  fs <- fs[!grepl(paste0(stem, "_Ai_"), fs)]
  if (!length(fs)) return(NULL)
  rbindlist(lapply(fs, read_gz), fill = TRUE) }
labels <- c("nbsr-hmc" = "NBSR (HMC)", "nbsr-laplace" = "NBSR (Laplace)", "nbsr-laplace-wideint" = "NBSR", "nbsr-laplace-eb" = "NBSR",
            "linda" = "LinDA", "aldex2" = "ALDEx2", "ancombc" = "ANCOM-BC", "ancombc2" = "ANCOM-BC2", "deseq2" = "DESeq2",
            "deseq2-shrink" = "DESeq2 + apeglm", "deseq2-poscounts" = "DESeq2 (poscounts)", "edger-qlf" = "edgeR QLF", "edger-lrt" = "edgeR LRT")
UNIT <- c("method", "simul_case", "sample_size", "replicate_id"); KEY <- c("simul_case", "sample_size", "replicate_id", "feature")

# Shared preparation: NBSR rows get the posterior-z Wald p and the lfsr; every method gets BH and Holm within a
# replicate over the features it tested; `call_own` is each method's stated rule.
prep <- function(D, ancombc2_rule = c("filter", "holm_filter")) {
  ancombc2_rule <- match.arg(ancombc2_rule)
  for (col in c("tested", "called", "passed_sens")) if (!col %in% names(D)) D[, (col) := NA]
  for (col in c("p_value", "q_value", "estimate", "se", "true_effect_abs", "p_dm")) if (!col %in% names(D)) D[, (col) := NA_real_]
  nb <- startsWith(D$method, "nbsr")
  D[nb, `:=`(tested = TRUE, estimate = post_mean, se = post_sd, p_value = 2 * pnorm(-abs(post_mean / post_sd)),
             lfsr = pmin(p_pos, 1 - p_pos))]
  D[!nb, lfsr := NA_real_]
  D[is.na(tested), tested := FALSE]; D[is.na(true_effect_abs), true_effect_abs := true_effect]
  D[tested == TRUE, q_bh := p.adjust(p_value, "BH"), by = UNIT]
  D[tested == TRUE, q_holm := p.adjust(p_value, "holm"), by = UNIT]
  D[tested == TRUE & !is.na(p_dm), q_dm := p.adjust(p_dm, "BH"), by = UNIT]
  D[, call_wald := tested & !is.na(q_bh) & q_bh <= ALPHA]
  D[, call_lfsr := nb & !is.na(lfsr) & lfsr <= LFSR_THR + LFSR_TOL]
  D[, call_dm := nb & !is.na(q_dm) & q_dm <= ALPHA]
  D[, call_own := fcase(nb, call_lfsr,
                        method == "ancombc2" & ancombc2_rule == "filter", tested & called %in% TRUE & passed_sens %in% TRUE,
                        method == "ancombc2", tested & !is.na(q_holm) & q_holm <= ALPHA & passed_sens %in% TRUE,
                        default = tested & called %in% TRUE)]
  D[is.na(call_own), call_own := FALSE]
  D[, stat := fifelse(tested & is.finite(estimate / se), abs(estimate / se), -Inf)]
  D }
confusion <- function(x, col, by) x[, .(TP = sum(get(col) & is_de), FP = sum(get(col) & !is_de), FN = sum(!get(col) & is_de),
                                       TN = sum(!get(col) & !is_de)), by = by][, `:=`(TPR = TP / (TP + FN), FDR = fifelse(TP + FP > 0, FP / (TP + FP), 0))][]
ap_one <- function(stat, de) { o <- order(-stat); de <- de[o]; tp <- cumsum(de); sum((tp / seq_along(de))[de]) / sum(de) }
power_at <- function(p, de, target) { o <- order(p); de <- de[o]; fdr <- cumsum(!de) / seq_along(de); ok <- which(fdr <= target)
  if (!length(ok)) return(0); sum(de[seq_len(max(ok))]) / sum(de) }
tables <- list()

## ---- Table 1: miRNA benchmark ------------------------------------------------------------------------------------
M <- prep(read_set("mirna"), "filter")
M1 <- M[method %in% c("nbsr-hmc", "nbsr-laplace", "linda", "aldex2", "ancombc", "ancombc2", "deseq2", "deseq2-shrink", "edger-qlf", "edger-lrt")]
t1 <- confusion(M1, "call_own", c("method", "sample_size"))
cov <- M1[tested == TRUE & is_de == TRUE & !is.na(q025), .(cov95_de = mean(q025 <= true_effect & true_effect <= q975)), by = .(method, sample_size)]
ap <- M1[, .(ap = ap_one(stat, is_de)), by = UNIT][, .(ap_median = median(ap), ap_mean = mean(ap)), by = .(method, sample_size)]
t1 <- Reduce(function(a, b) merge(a, b, by = c("method", "sample_size"), all.x = TRUE), list(t1, cov, ap))
t1[startsWith(method, "edger"), cov95_de := NA]      # edgeR provides no interval
t1[, label := labels[method]]; setorder(t1, sample_size, method)
tables$table1_mirna <- t1[, .(label, method, sample_size, TPR = round(TPR, 3), FDR = round(FDR, 3), cov95_de = round(cov95_de, 3), AP_median = round(ap_median, 3), TP, FP, FN)]

## ---- Table 2 and the appendix grid: microbiome benchmark ---------------------------------------------------------
B <- prep(read_set("microbiome"), "holm_filter")
four <- c("nbsr-laplace-wideint", "linda", "aldex2", "ancombc2")
B4 <- B[method %in% four]
common <- B4[, .(ok = all(tested) && uniqueN(method) == 4), by = KEY][ok == TRUE, ..KEY]
C <- merge(B4, common, by = KEY)
C[, call_rec := fifelse(startsWith(method, "nbsr"), call_wald, call_own)]
t2 <- rbind(confusion(C, "call_rec", c("method", "simul_case", "sample_size"))[, rule := "stated rule"],
            confusion(C[method == "nbsr-laplace-wideint"], "call_lfsr", c("method", "simul_case", "sample_size"))[, rule := "NBSR, lfsr <= 1e-3"],
            confusion(C[method == "nbsr-laplace-wideint" & !is.na(p_dm)], "call_dm", c("method", "simul_case", "sample_size"))[, rule := "NBSR, delta-method shift + BH"])
t2[, label := labels[method]]; setorder(t2, simul_case, sample_size, rule, method)
tables$table2_microbiome <- t2[, .(rule, label, simul_case, sample_size, power = round(TPR, 3), FDR = round(FDR, 3), TP, FP, FN)]
mf <- rbind(C[, .(statistic = "own p-value", tpr05 = power_at(fifelse(tested, p_value, 1), is_de, 0.05),
                  tpr10 = power_at(fifelse(tested, p_value, 1), is_de, 0.10)), by = .(method, simul_case, sample_size)],
            C[method == "nbsr-laplace-wideint", .(statistic = "lfsr", tpr05 = power_at(lfsr, is_de, 0.05), tpr10 = power_at(lfsr, is_de, 0.10)), by = .(method, simul_case, sample_size)],
            C[method == "nbsr-laplace-wideint" & !is.na(p_dm), .(statistic = "delta-method", tpr05 = power_at(p_dm, is_de, 0.05), tpr10 = power_at(p_dm, is_de, 0.10)), by = .(method, simul_case, sample_size)])
mf[, label := labels[method]]; setorder(mf, simul_case, sample_size, method, statistic)
tables$matched_fdr <- mf[, .(label, statistic, simul_case, sample_size, power_at_fdr05 = round(tpr05, 3), power_at_fdr10 = round(tpr10, 3))]
# calibration and fold-change accuracy of the four methods on the common set
cal <- C[tested == TRUE & !is.na(q025), .(bias = mean(estimate - true_effect), rmse = sqrt(mean((estimate - true_effect)^2)),
                                          cov95 = mean(q025 <= true_effect & true_effect <= q975)),
         by = .(method, simul_case, sample_size, truth = fifelse(is_de, "DE", "null"))]
cal[, label := labels[method]]; setorder(cal, simul_case, sample_size, method, truth)
tables$microbiome_calibration <- cal[, .(label, simul_case, sample_size, truth, bias = round(bias, 3), rmse = round(rmse, 3), cov95 = round(cov95, 3))]

## ---- per-sample dispersion variants (fixed and drawn concentration) --------------------------------------------
Ai <- read_set("microbiome_Ai")
disp <- rbind(if (!is.null(Ai)) prep(Ai)[, generator := "kappa_i drawn"], B[method %in% c("nbsr-laplace-wideint", "nbsr-laplace-wsample_wideint", "nbsr-laplace-wdepth_wideint")][, generator := "kappa fixed"], fill = TRUE)
if (nrow(disp)) {
  dv <- confusion(disp, "call_wald", c("generator", "method", "simul_case", "sample_size"))
  extra <- disp[, .(null_bias = mean(post_mean[!is_de]), cov95_de = mean((q025 <= true_effect & true_effect <= q975)[is_de])), by = .(generator, method, simul_case, sample_size)]
  dv <- merge(dv, extra, by = c("generator", "method", "simul_case", "sample_size")); setorder(dv, generator, simul_case, method, sample_size)
  tables$dispersion_variants <- dv[, .(generator, method, simul_case, sample_size, power = round(TPR, 3), FDR = round(FDR, 3), null_bias = round(null_bias, 3), cov95_de = round(cov95_de, 3))] }

## ---- Table 3: SELEX ----------------------------------------------------------------------------------------------
X <- prep(read_set("selex"), "holm_filter")
X[, up := !is.na(estimate) & estimate > 0]
sx <- function(col, rule) X[, .(on27 = sum(get(col) & nixon27), others_up = sum(get(col) & up & !nixon27), others_down = sum(get(col) & !up & !nixon27)),
                             by = UNIT][, lapply(.SD, mean), by = .(method, sample_size), .SDcols = c("on27", "others_up", "others_down")][, rule := rule]
t3 <- rbind(sx("call_own", "own rule")[!startsWith(method, "nbsr") | TRUE], X[startsWith(method, "nbsr")][, .(on27 = sum(call_wald & nixon27), others_up = sum(call_wald & up & !nixon27), others_down = sum(call_wald & !up & !nixon27)), by = UNIT][, lapply(.SD, mean), by = .(method, sample_size), .SDcols = c("on27", "others_up", "others_down")][, rule := "NBSR, posterior-z Wald + BH"],
            X[startsWith(method, "nbsr") & !is.na(p_dm)][, .(on27 = sum(call_dm & nixon27), others_up = sum(call_dm & up & !nixon27), others_down = sum(call_dm & !up & !nixon27)), by = UNIT][, lapply(.SD, mean), by = .(method, sample_size), .SDcols = c("on27", "others_up", "others_down")][, rule := "NBSR, delta-method shift + BH"])
t3[startsWith(method, "nbsr") & rule == "own rule", rule := "NBSR, lfsr <= 1e-3"]
setorder(t3, sample_size, method, rule)
tables$table3_selex <- t3[, .(rule, method, sample_size, on27 = round(on27, 1), others_up = round(others_up, 1), others_down = round(others_down, 1))]

## ---- the label-permutation null ------------------------------------------------------------------------------------
P <- prep(read_set("permutation"), "holm_filter")
P[, call := fifelse(startsWith(method, "nbsr"), call_wald, call_own)]
pr <- P[, .(tested = sum(tested), fp = sum(call), fp_lfsr = sum(call_lfsr), p05 = mean(p_value[tested] < 0.05, na.rm = TRUE)), by = UNIT]
perm <- pr[, .(replicates = .N, fp_mean = mean(fp), pct_with_call = 100 * mean(fp > 0), worst = max(fp), fpr = sum(fp) / sum(tested), fpr_lfsr = sum(fp_lfsr) / sum(tested)), by = .(method, simul_case, sample_size)]
raw <- pr[, .(raw_p_below_05 = mean(p05, na.rm = TRUE)), by = method]
perm <- merge(perm, raw, by = "method"); perm[, label := labels[method]]; setorder(perm, simul_case, sample_size, method)
tables$permutation <- perm[, .(label, simul_case, sample_size, fp_mean = round(fp_mean, 2), pct_with_call = round(pct_with_call), worst, fpr = signif(fpr, 2), fpr_lfsr = signif(fpr_lfsr, 2), raw_p_below_05 = round(raw_p_below_05, 3))]

## ---- the studies with their real labels ----------------------------------------------------------------------------
RL <- prep(read_set("real_labels"), "holm_filter")
RL[, call := fifelse(startsWith(method, "nbsr"), call_wald, call_own)]
EST <- c("Fusobacterium", "Parvimonas", "Porphyromonas", "Peptostreptococcus", "Gemella", "Solobacterium")
RL[, established := simul_case == "BAXTER_CRC" & feature %in% EST]
rl <- RL[, .(tested = sum(tested), up = sum(call & estimate > 0), down = sum(call & estimate < 0), up_lfsr = sum(call_lfsr & estimate > 0), down_lfsr = sum(call_lfsr & estimate < 0),
             established_up = sum(call & established & estimate > 0), established_tested = sum(tested & established)), by = .(method, simul_case)]
rl[, label := labels[method]]; setorder(rl, simul_case, method)
tables$real_labels <- rl[, .(label, simul_case, tested, up, down, up_lfsr, down_lfsr, established_up, established_tested)]
tables$established_genera <- RL[established == TRUE, .(method, genus = feature, tested, estimate = round(estimate, 2), q = signif(fifelse(startsWith(method, "nbsr"), q_bh, fifelse(method == "ancombc2", q_holm, q_value)), 2), called = call)][order(genus, method)]

## ---- appendix: calibration, decision rules, Laplace against HMC, Wald rules ----------------------------------------
H <- M[method %in% c("nbsr-hmc", "nbsr-laplace", "nbsr-laplace-eb")]
H[, `:=`(err = post_mean - true_effect, z = (post_mean - true_effect) / post_sd, w95 = q975 - q025,
         c50 = q25 <= true_effect & true_effect <= q75, c90 = q05 <= true_effect & true_effect <= q95,
         c95 = q025 <= true_effect & true_effect <= q975, c99 = q005 <= true_effect & true_effect <= q995,
         call95 = q025 > 0 | q975 < 0, ppmax = pmax(p_pos, 1 - p_pos))]
calib <- function(x) x[, .(n = .N, cov50 = mean(c50), cov90 = mean(c90), cov95 = mean(c95), cov99 = mean(c99), z_sd = sd(z), rmse = sqrt(mean(err^2)), bias = mean(err))]
cb <- rbind(H[, c(list(features = "all"), calib(.SD)), by = .(method, sample_size)], H[, c(list(features = fifelse(is_de[1], "DE", "null")), calib(.SD)), by = .(method, sample_size, is_de)][, is_de := NULL])
setorder(cb, method, sample_size, features)
tables$calibration <- cb[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]
bfdr <- function(lfsr, alpha) { o <- order(lfsr); k <- max(c(0L, which(cumsum(lfsr[o]) / seq_along(o) <= alpha))); call <- logical(length(lfsr)); if (k > 0) call[o[seq_len(k)]] <- TRUE; call }
H[, bfdr01 := bfdr(lfsr, 0.01), by = UNIT]
rules <- list("95% interval excludes 0" = "call95", "lfsr <= 0.01" = "r01", "lfsr <= 0.001" = "call_lfsr", "mean lfsr <= 0.01 (Bayesian FDR)" = "bfdr01")
H[, r01 := ppmax >= 0.99]
dr <- rbindlist(lapply(names(rules), function(r) confusion(H, rules[[r]], c("method", "sample_size"))[, rule := r]))
dr[, FPR := FP / (FP + TN)]; setorder(dr, method, sample_size)
tables$decision_rules <- dr[, .(method, sample_size, rule, TP, FP, FN, TPR = round(TPR, 3), FPR = round(FPR, 4), FDR = round(FDR, 3))]
sep <- H[, .(perfect = min(ppmax[is_de]) > max(ppmax[!is_de]), auc = { s1 <- ppmax[is_de]; s0 <- ppmax[!is_de]; mean(outer(s1, s0, ">") + 0.5 * outer(s1, s0, "==")) }), by = UNIT][
  , .(perfect_separation = sum(perfect), auc_min = min(auc), auc_median = median(auc)), by = .(method, sample_size)]
lv <- merge(dcast(cb[features != "all", .(method, sample_size, features, cov95, z_sd)], method + sample_size ~ features, value.var = c("cov95", "z_sd")),
            merge(dr[rule == "lfsr <= 0.001", .(method, sample_size, TP, FP, FDR)], sep, by = c("method", "sample_size")), by = c("method", "sample_size"))
tables$laplace_vs_hmc <- lv[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]
W <- H[method %in% c("nbsr-hmc", "nbsr-laplace") & !is.na(p_dm)]
wr <- rbind(confusion(W, "call_dm", c("method", "sample_size"))[, rule := "delta-method Wald + BH"], confusion(W, "call_wald", c("method", "sample_size"))[, rule := "posterior-z Wald + BH"],
            confusion(W, "call_lfsr", c("method", "sample_size"))[, rule := "lfsr <= 0.001"])
setorder(wr, method, sample_size, rule); tables$wald_rules <- wr[, .(method, sample_size, rule, TP, FP, TPR = round(TPR, 3), FDR = round(FDR, 3))]

## ---- appendix: the ablation variants --------------------------------------------------------------------------------
AB <- read_set("ablation")
if (!is.null(AB)) {
  AB <- prep(AB); AB[, `:=`(z = (post_mean - true_effect) / post_sd, c95 = q025 <= true_effect & true_effect <= q975, ppmax = pmax(p_pos, 1 - p_pos))]
  base <- H[method == "nbsr-hmc"][, c("method", names(AB)[names(AB) %in% names(H)]), with = FALSE]
  AB <- rbind(AB, H[method == "nbsr-hmc", intersect(names(AB), names(H)), with = FALSE], fill = TRUE)
  ab <- merge(AB[, .(units = uniqueN(replicate_id), cov95_de = mean(c95[is_de]), zsd_de = sd(z[is_de]), bias_de = mean((post_mean - true_effect)[is_de]), rmse_de = sqrt(mean(((post_mean - true_effect)[is_de])^2)),
                     cov95_null = mean(c95[!is_de]), bias_null = mean(post_mean[!is_de]), width95_null = median((q975 - q025)[!is_de]),
                     TP = sum(call_lfsr & is_de), FP = sum(call_lfsr & !is_de), TPR = sum(call_lfsr & is_de) / sum(is_de), FDR = sum(call_lfsr & !is_de) / max(1, sum(call_lfsr))), by = .(method, sample_size)],
              AB[, .(perfect = min(ppmax[is_de]) > max(ppmax[!is_de]), auc = { s1 <- ppmax[is_de]; s0 <- ppmax[!is_de]; mean(outer(s1, s0, ">") + 0.5 * outer(s1, s0, "==")) }), by = UNIT][
                , .(perfect_separation = sum(perfect), auc_min = min(auc)), by = .(method, sample_size)], by = c("method", "sample_size"))
  setorder(ab, sample_size, method); tables$ablations <- ab[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)] }

## ---- write ----------------------------------------------------------------------------------------------------------
md <- function(x) { x <- as.data.frame(x); c(paste0("| ", paste(names(x), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(x)), collapse = "|"), "|"), apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) }
lines <- c("# Tables rebuilt from the released outputs", "", paste("data:", data_dir), "")
for (nm in names(tables)) { fwrite(tables[[nm]], file.path(out, paste0(nm, ".csv"))); lines <- c(lines, paste("##", nm), "", md(tables[[nm]]), "") }
writeLines(lines, file.path(out, "tables.md"))
cat("wrote", length(tables), "tables to", out, "\n")
