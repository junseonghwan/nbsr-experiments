# Diagnostic plots for the Dirichlet-multinomial fit to Baxter healthy samples.
suppressPackageStartupMessages(library(data.table))
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
set.seed(1)
args <- commandArgs(trailingOnly = TRUE)
ds <- if (length(args) >= 1) args[1] else "crc_baxter_genus"
ctrl <- if (length(args) >= 2) args[2] else "H"
d <- file.path("data/microbiome/processed", ds)
Y <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE))
X <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
est <- fread(file.path(d, paste0("dirichlet_", ctrl, ".csv")))
cts <- Y[est$feature, X$trt == ctrl]
a <- est$alpha_dm_mle; A <- sum(a); m <- a / A
n <- colSums(cts); N <- ncol(cts); K <- nrow(cts)
p_obs <- sweep(cts, 2, n, "/")

# ---- simulate replicate datasets from the fitted DM with the observed library sizes
sim_dm <- function() {
  g <- matrix(rgamma(K * N, shape = a, rate = 1), K, N)          # Dirichlet draws by column
  p <- sweep(g, 2, colSums(g), "/")
  sapply(seq_len(N), function(i) rmultinom(1, n[i], p[, i]))
}
R <- 200
sims <- replicate(R, sim_dm(), simplify = FALSE)

# per-genus statistics: observed vs simulated (median and 95% band)
stat_genus <- function(M) cbind(mean = rowMeans(sweep(M, 2, colSums(M), "/")),
                                prev = rowMeans(M > 0),
                                var  = apply(sweep(M, 2, colSums(M), "/"), 1, var))
obs_g <- stat_genus(cts)
sim_g <- simplify2array(lapply(sims, stat_genus))       # K x 3 x R
q_g <- function(j, p) apply(sim_g[, j, ], 1, quantile, p)

# per-sample statistics
stat_sample <- function(M) { P <- sweep(M, 2, colSums(M), "/")
  cbind(richness = colSums(M > 0), shannon = -colSums(ifelse(P > 0, P * log(P), 0)), top1 = apply(P, 2, max)) }
obs_s <- stat_sample(cts); sim_s <- simplify2array(lapply(sims, stat_sample))

png(file.path(d, paste0("dirichlet_diagnostics_", ctrl, ".png")), width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))
# 1 mean proportion
plot(m, obs_g[, "mean"], log = "xy", pch = 16, col = rgb(0,0,0,.5), xlab = "fitted alpha_j / A", ylab = "observed mean proportion",
     main = "Mean proportion"); abline(0, 1, col = 2)
# 2 prevalence with simulated band
o <- order(m); plot(m[o], obs_g[o, "prev"], log = "x", pch = 16, col = rgb(0,0,0,.5), ylim = c(0, 1),
     xlab = "fitted alpha_j / A", ylab = "fraction of samples with count > 0", main = "Prevalence (band: 95% of DM replicates)")
polygon(c(m[o], rev(m[o])), c(q_g("prev", .025)[o], rev(q_g("prev", .975)[o])), col = rgb(1,0,0,.2), border = NA)
lines(m[o], q_g("prev", .5)[o], col = 2)
# 3 variance of proportions
plot(q_g("var", .5), obs_g[, "var"], log = "xy", pch = 16, col = rgb(0,0,0,.5),
     xlab = "DM-replicate median variance of proportion", ylab = "observed variance of proportion", main = "Between-sample variance")
abline(0, 1, col = 2); abline(log10(1/3), 1, col = 2, lty = 3); legend("topleft", c("identity", "1/3 of predicted"), lty = c(1, 3), col = 2, bty = "n")
# 4 richness
hist(sim_s[, 1, ], breaks = 40, col = rgb(1,0,0,.3), border = NA, freq = FALSE,
     main = "Per-sample richness (genera > 0)", xlab = "richness"); lines(density(obs_s[, "richness"]), lwd = 2)
legend("topright", c("DM replicates", "observed"), fill = c(rgb(1,0,0,.3), NA), border = NA, lty = c(NA, 1), lwd = c(NA, 2), bty = "n")
# 5 shannon
hist(sim_s[, 2, ], breaks = 40, col = rgb(1,0,0,.3), border = NA, freq = FALSE, main = "Per-sample Shannon diversity", xlab = "Shannon")
lines(density(obs_s[, "shannon"]), lwd = 2)
# 6 dominance
hist(sim_s[, 3, ], breaks = 40, col = rgb(1,0,0,.3), border = NA, freq = FALSE, main = "Per-sample share of the top genus", xlab = "max proportion")
lines(density(obs_s[, "top1"]), lwd = 2)
dev.off()

cat(sprintf("Per-sample summaries, observed median [IQR] vs DM-replicate median [IQR]:\n"))
for (j in 1:3) cat(sprintf("  %-9s obs %.2f [%.2f, %.2f]   sim %.2f [%.2f, %.2f]\n", colnames(obs_s)[j],
  median(obs_s[, j]), quantile(obs_s[, j], .25), quantile(obs_s[, j], .75),
  median(sim_s[, j, ]), quantile(sim_s[, j, ], .25), quantile(sim_s[, j, ], .75)))
cat(sprintf("Genera whose observed prevalence falls outside the 95%% DM band: %d of %d\n",
            sum(obs_g[, "prev"] < q_g("prev", .025) | obs_g[, "prev"] > q_g("prev", .975)), K))
cat(sprintf("Genera whose observed variance falls outside the 95%% DM band: %d of %d (below: %d)\n",
            sum(obs_g[, "var"] < q_g("var", .025) | obs_g[, "var"] > q_g("var", .975)), K, sum(obs_g[, "var"] < q_g("var", .025))))
cat("Plot:", file.path(d, paste0("dirichlet_diagnostics_", ctrl, ".png")), "\n")
