# Estimate the lognormal layer of the logistic-normal Dirichlet-multinomial generator (simulate_lndm.R)
# from a template's control samples.
#
# Usage: Rscript data_processor/estimate_lognormal.R [--template crc_baxter_genus] [--control H] [--rank 5] [--A 100]
#          [--diag hybrid|moment|pln] [--cap 3]
#
# VALIDATED SETTINGS for crc_baxter_genus (2026-10-08, data_processor/validate_generator.R against the 172 healthy
# samples): --A 100 --rank 5 --diag hybrid --cap 3. Correlation spread, strong pairs and leading components match
# the template within 0.02; prevalence and the MED/LOW zero fractions match; the median per-genus dispersion ratio
# matches. Known shortfall: the extremes of per-genus dispersion are compressed (10th/90th percentile of the
# variance ratio 0.36/1.9 against 0.02/3.9), because a single Dirichlet concentration bounds how stable a genus can
# be and the cap bounds how bursty. The DM maximum-likelihood A (38) reproduces zeros but leaves no room for the
# lognormal diagonal; A >= 300 loses the sparse genera's zeros.
#
# Model the generator draws from, per sample i:
#   Z_i ~ N(0, Sigma),  alpha_i = A * softmax(log alpha + Z_i),  p_i ~ Dirichlet(alpha_i),  y_i ~ Multinomial(n_i, p_i)
# with alpha the Dirichlet-multinomial fit (dirichlet_<control>.csv) and
#   Sigma = L L' + diag(d):  a rank-q factor part estimated among the PREVALENT genera (present in >= 50% of
#   control samples; rare genera get zero loadings because their correlation is not estimable), and a genus-
#   specific diagonal d_j. Three estimators of the diagonal (hybrid is the default):
#     hybrid            moment matching for the prevalent genera, whose variance is measurable, and the PLN
#                       remainder for the rest, whose only measurable property is how often they are zero and
#                       which the PLN variance reproduces.
#     moment (default)  match each prevalent genus's observed variance of proportion in the template: with total
#                       log-scale variance s_j^2 = (L L')_jj + d_j the lognormal layer adds ~ pi_j^2 (exp(s_j^2) - 1)
#                       to the Dirichlet-multinomial variance, so s_j^2 = log(1 + (v_obs - v_dm(A)) / pi_j^2),
#                       floored at 0; rare genera (not estimable) get d_j = 0.
#     pln               d_j = max(0, tau2_j - trigamma(A pi_j) - (L L')_jj) from the PLN latent variance; this
#                       over-disperses sparse genera, whose PLN variance is dominated by zeros.
#   The simulator centres the layer (mean log alpha_j - s_j^2 / 2) so each genus keeps its template mean share.
#   Sigma = 0 recovers the plain Dirichlet-multinomial exactly.
#
# The Poisson-lognormal (PLNmodels) is fitted to all genera with the library size as offset; its full covariance is
# used for tau2 and, restricted to the prevalent genera, for the factor part. A is a knob: the DM maximum-likelihood
# value by default (continuity with the DM benchmark); larger A moves variance from the Dirichlet layer into d.
#
# Output: <processed>/<template>/lognormal_<control>_A<A>_q<rank>.rds with alpha, A, L, d, q, genera, prevalent,
# tau2, and the PLN loglik, plus a short diagnostics print.
suppressPackageStartupMessages({library(data.table); library(PLNmodels)})
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
opt <- function(flag, default) { a <- commandArgs(TRUE); i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
template <- opt("--template", "crc_baxter_genus"); control <- opt("--control", "H")
q <- as.integer(opt("--rank", "5")); A_opt <- opt("--A", "100"); diag_mode <- opt("--diag", "hybrid")
cap <- as.numeric(opt("--cap", "3"))   # a lognormal variance above ~3 lets single genera take over samples ("blooms")
stopifnot(diag_mode %in% c("hybrid", "moment", "pln"))
set.seed(1)
d <- file.path("data/microbiome/processed", template)
est <- fread(file.path(d, paste0("dirichlet_", control, ".csv")))
Y <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE)); X <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
H <- t(Y[est$feature, X$trt == control]); lib <- rowSums(H)
alpha <- setNames(est$alpha_dm_mle, est$feature); A_dm <- sum(alpha); A <- if (A_opt == "dm") A_dm else as.numeric(A_opt); pi0 <- alpha / A_dm
cat(sprintf("Template %s (%s): %d samples x %d genera; A_DM = %.1f, A used = %.1f, rank q = %d\n", template, control, nrow(H), ncol(H), A_dm, A, q))
fit <- PLN(H ~ 1 + offset(log(lib)), control = PLN_param(trace = 0))
Sig <- sigma(fit); tau2 <- diag(Sig); genera <- colnames(H)
prevalent <- colMeans(H > 0) >= 0.5
# factor part among prevalent genera: top-q eigenpairs of their PLN covariance, rare genera get zero loadings
Sp <- Sig[prevalent, prevalent]; e <- eigen(Sp, symmetric = TRUE)
L <- matrix(0, length(genera), q, dimnames = list(genera, NULL)); L[prevalent, ] <- e$vectors[, 1:q] %*% diag(sqrt(pmax(e$values[1:q], 0)), q)
# diagonal: PLN latent variance minus the Dirichlet layer's own log-scale variance, minus what the factors already carry
dir_var <- trigamma(A * pi0); fac_var <- rowSums(L^2)
d_pln <- pmax(0, tau2 - dir_var - fac_var)
Pm <- H / lib; v_obs <- apply(Pm, 2, var)
v_dm <- pi0 * (1 - pi0) / (A + 1) + pi0 * (1 - pi0) * mean(1 / lib)
s2 <- ifelse(v_obs > v_dm, log(1 + (v_obs - v_dm) / pi0^2), 0)
d_mom <- ifelse(prevalent, pmax(0, s2 - fac_var), 0)
dvec <- switch(diag_mode, pln = d_pln, moment = d_mom, hybrid = ifelse(prevalent, d_mom, d_pln))
dvec <- pmax(0, pmin(dvec, cap - fac_var))      # cap the total s2 = factor + diagonal
s2_total <- fac_var + dvec
cat(sprintf("factor part: %d prevalent genera, top-%d eigenvalues carry %.0f%% of their PLN covariance trace\n", sum(prevalent), q, 100 * sum(e$values[1:q]) / sum(e$values)))
strat <- fifelse(est$mean_prop > 0.01, "HIGH", fifelse(est$mean_prop > 0.001, "MED", "LOW"))
sm <- data.table(stratum = strat, tau2 = tau2, dirichlet = dir_var, factor = fac_var, diag = dvec, prevalent = prevalent)
cat(sprintf("diagonal estimator: %s\n", diag_mode))
cat("\nvariance split by stratum (medians): PLN latent tau2 | Dirichlet layer at this A | factor part | lognormal diagonal\n")
print(sm[, .(n = .N, tau2 = round(median(tau2), 2), dirichlet = round(median(dirichlet), 2), factor = round(median(factor), 2), diag = round(median(diag), 2), share_with_diag_gt0 = round(mean(diag > 0), 2)), by = stratum])
out <- file.path(d, sprintf("lognormal_%s_A%.0f_q%d_%s%s.rds", control, A, q, diag_mode, if (is.finite(cap)) sprintf("_cap%g", cap) else ""))
saveRDS(list(template = template, control = control, alpha = alpha, A = A, A_dm = A_dm, q = q, L = L, d = dvec, s2 = s2_total,
             diag_mode = diag_mode, genera = genera, prevalent = prevalent, tau2 = tau2, pln_loglik = fit$loglik, Sigma_pln = Sig), out)
cat("\nsaved", out, "\n")
