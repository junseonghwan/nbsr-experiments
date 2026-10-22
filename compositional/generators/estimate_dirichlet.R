# Estimate Dirichlet / Dirichlet-multinomial concentration parameters for the
# the control group of a processed dataset. Usage: Rscript estimate_dirichlet.R <processed dir> <control label>
# Three estimators are compared:
#   (1) estimate_dirichlet():        Dirichlet MLE on observed proportions (miRNA paper approach, 1e-6 smoothing)
#   (2) optimize_concentration():    two-step DM: alpha_bar = mean proportions, precision s by fixed point
#   (3) DM MLE:                      full Dirichlet-multinomial likelihood on counts, optim over log alpha
suppressPackageStartupMessages(library(data.table))
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source("data_processor/dirichlet_estimation_functions.R")   # copied from nbsr-experiments/paper/ (published miRNA paper)
normalize <- function(x) x / sum(x)

args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args) >= 1) args[1] else "crc_baxter_genus"
ctrl <- if (length(args) >= 2) args[2] else "H"
d <- file.path("data/microbiome/processed", ds)
Y <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE))   # features x samples
X <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(Y), rownames(X)))
cts <- Y[, X$trt == ctrl]
cat("Control group", ctrl, ":", ncol(cts), "samples,", nrow(cts), "genera; library size median", median(colSums(cts)), "\n")
keep <- rowSums(cts) > 0
cat("Genera with zero counts in the healthy group (dropped):", sum(!keep), "\n")
cts <- cts[keep, ]
K <- nrow(cts); N <- ncol(cts)

# (1) Dirichlet on proportions, as in estimate_microRNAome.R
t0 <- Sys.time()
ret1 <- estimate_dirichlet(t(cts), rep(1, K), max_iters = 1000)
a1 <- ret1$alpha_hat; names(a1) <- rownames(cts)
cat(sprintf("(1) Dirichlet-on-proportions: converged=%s, A=%.2f, %.1fs\n", ret1$converged, sum(a1), as.numeric(Sys.time()-t0, units="secs")))

# (2) two-step DM
alpha_bar <- rowMeans(apply(cts, 2, normalize))
ret2 <- optimize_concentration(cts, alpha_bar, s_init = 1, max_iter = 10000)
a2 <- alpha_bar * ret2$s
cat(sprintf("(2) Two-step DM (mean proportions x precision): s=A=%.2f after %d iterations\n", ret2$s, ret2$iters))

# (3) full DM MLE on counts
dm_nll <- function(la) {
  a <- exp(la); A <- sum(a); n <- colSums(cts)
  -(N * lgamma(A) - sum(lgamma(n + A)) + sum(lgamma(cts + a) - lgamma(a)))
}
dm_grad <- function(la) {
  a <- exp(la); A <- sum(a); n <- colSums(cts)
  g <- N * digamma(A) - sum(digamma(n + A)) + rowSums(digamma(cts + a)) - N * digamma(a)
  -g * a
}
t0 <- Sys.time()
opt <- optim(log(a2), dm_nll, dm_grad, method = "L-BFGS-B", control = list(maxit = 20000, factr = 1e2, pgtol = 1e-8))
a3 <- exp(opt$par); names(a3) <- rownames(cts)
cat(sprintf("(3) DM MLE: convergence code %d, A=%.2f, negloglik=%.1f, %.1fs\n", opt$convergence, sum(a3), opt$value, as.numeric(Sys.time()-t0, units="secs")))
cat(sprintf("    negloglik at two-step estimate: %.1f; at Dirichlet-on-proportions estimate: %.1f\n", dm_nll(log(a2)), dm_nll(log(pmax(a1, 1e-12)))))

# predictive checks for the DM MLE: zero fraction and proportion variance
n <- colSums(cts); A <- sum(a3)
pzero <- sapply(seq_len(N), function(i) exp(lgamma(n[i] + A - a3) - lgamma(n[i] + A) + lgamma(A) - lgamma(A - a3)))
cat(sprintf("    observed zero fraction %.3f vs DM-predicted %.3f\n", mean(cts == 0), mean(pzero)))
p_obs <- apply(cts, 2, normalize)
v_obs <- apply(p_obs, 1, var); m <- a3 / A
v_dm <- m * (1 - m) / (A + 1) + mean(1 / n) * m * (1 - m) * A / (A + 1)   # Dirichlet + multinomial sampling
cat(sprintf("    proportion variance, observed vs DM-predicted (log10 ratio median over genera with mean prop > 1e-3): %.2f\n",
            median(log10(v_obs[m > 1e-3] / v_dm[m > 1e-3]))))

out <- data.table(feature = rownames(cts), mean_prop = m, prevalence = rowMeans(cts > 0),
                  alpha_dirichlet_props = a1, alpha_dm_twostep = a2, alpha_dm_mle = a3)[order(-mean_prop)]
fwrite(out, file.path(d, paste0("dirichlet_", ctrl, ".csv")))
cat("\nAlpha (DM MLE) quantiles:\n"); print(round(quantile(a3, c(0, .1, .25, .5, .75, .9, 1)), 4))
cat("count alpha > 1:", sum(a3 > 1), "; alpha in (0.1,1]:", sum(a3 > .1 & a3 <= 1), "; alpha <= 0.1:", sum(a3 <= .1), "\n")
cat("\nTop 10 genera:\n"); print(out[1:10, .(feature, mean_prop = round(mean_prop, 4), prevalence = round(prevalence, 2), alpha_dm_mle = round(alpha_dm_mle, 3), alpha_dirichlet_props = round(alpha_dirichlet_props, 3))])
cat("\nWritten:", file.path(d, paste0("dirichlet_", ctrl, ".csv")), "\n")
