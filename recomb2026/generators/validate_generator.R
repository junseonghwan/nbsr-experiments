# Does a simulated benchmark look like the template it was anchored on? Four checks on the control samples:
# between-genus correlation (spread, strong pairs, leading components), per-genus dispersion relative to a
# single-concentration Dirichlet-multinomial, zero fraction by stratum, and mean proportions.
#
# Usage: Rscript data_processor/validate_generator.R <simulated root> [--template crc_baxter_genus] [--control H]
#          [--case ALL25] [--size 80] [--nsamples 172]
# Compares the template's control samples with <nsamples> control samples drawn from the simulated replicates of
# one case and size (controls are unperturbed, so any case serves). Prints a table: template vs simulated.
suppressPackageStartupMessages(library(data.table)); set.seed(1)
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
a <- commandArgs(TRUE); sim_root <- a[1]
opt <- function(flag, default) { i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
template <- opt("--template", "crc_baxter_genus"); control <- opt("--control", "H"); case <- opt("--case", "ALL25")
size <- opt("--size", "80"); nsamp <- as.integer(opt("--nsamples", "172"))
d <- file.path("data/microbiome/processed", template)
est <- fread(file.path(d, paste0("dirichlet_", control, ".csv"))); alpha <- setNames(est$alpha_dm_mle, est$feature); A <- sum(alpha)
Y <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE)); X <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
real <- Y[est$feature, X$trt == control]
# simulated controls: pool the "null" samples over replicates, keep the template's genera (absent ones are zero rows)
reps <- Sys.glob(file.path(sim_root, case, paste0("sample_", size), "rep*"))
sim <- do.call(cbind, lapply(reps, function(r) { Ys <- as.matrix(read.csv(file.path(r, "Y.csv"), row.names = 1, check.names = FALSE)); Xs <- read.csv(file.path(r, "X.csv"), row.names = 1, check.names = FALSE)
  Ys <- Ys[, Xs$trt == "null", drop = FALSE]; M <- matrix(0, length(alpha), ncol(Ys), dimnames = list(names(alpha), paste0(basename(r), "_", colnames(Ys)))); M[rownames(Ys), ] <- Ys; M }))
sim <- sim[, sample(ncol(sim), min(nsamp, ncol(sim)))]
strat <- fifelse(est$mean_prop > 0.01, "HIGH", fifelse(est$mean_prop > 0.001, "MED", "LOW")); names(strat) <- est$feature
clr_cor <- function(M, keep) { P <- sweep(M[keep, , drop = FALSE] + 0.5, 2, colSums(M[keep, , drop = FALSE] + 0.5), "/"); L <- log(P); C <- L - rep(colMeans(L), each = nrow(L)); cor(t(C)) }
diag_one <- function(M, label) {
  keep <- rowMeans(M > 0) >= 0.5; R <- clr_cor(M, keep); r <- R[upper.tri(R)]; ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  P <- sweep(M, 2, colSums(M), "/"); pi0 <- alpha / A; v <- apply(P, 1, var); vdm <- pi0 * (1 - pi0) / (A + 1) + pi0 * (1 - pi0) * mean(1 / colSums(M)); ratio <- (v / vdm)[keep]
  zf <- tapply(rowMeans(M == 0), strat[rownames(M)], mean)
  data.table(dataset = label, samples = ncol(M), prevalent_units = sum(keep), corr_sd = sd(r), pairs_r_gt_0.5 = sum(abs(r) > 0.5), pc1 = ev[1] / sum(ev), pc5 = sum(ev[1:5]) / sum(ev),
             var_ratio_q10 = quantile(ratio, .1), var_ratio_med = median(ratio), var_ratio_q90 = quantile(ratio, .9),
             zero_HIGH = zf[["HIGH"]], zero_MED = zf[["MED"]], zero_LOW = zf[["LOW"]], lib_median = median(colSums(M)))
}
out <- rbind(diag_one(real, "template"), diag_one(sim, basename(sim_root)))
cat(sprintf("template %s (%s) vs simulated %s [%s/sample_%s, %d control samples]\n", template, control, sim_root, case, size, ncol(sim)))
print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
cat("\nrows: correlation spread and strong pairs among prevalent genera (CLR, pseudocount 0.5); PC shares of the correlation matrix;\n")
cat("per-genus variance of proportion over the single-A Dirichlet prediction (quantiles over prevalent genera; DM gives ~1 everywhere);\n")
cat("zero fraction by stratum; median library size.\n")
