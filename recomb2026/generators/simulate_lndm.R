# Logistic-normal Dirichlet-multinomial simulator, anchored on a fitted template (default: Baxter 2016 healthy
# stool, genus level). Per sample i,
#   Z_i ~ N(0, Sigma),  alpha_i = A * softmax(log alpha - diag(Sigma)/2 + Z_i)  [or A_i = sum(...) with --vary_A 1],
#   p_i ~ Dirichlet(alpha_i),  y_i ~ Multinomial(n_i, p_i).
# The -diag(Sigma)/2 term is the lognormal mean correction: each genus keeps the template's mean share in
# expectation, so the layer changes the spread and not the composition.
# Sigma = L L' + diag(d) comes from data_processor/estimate_lognormal.R. The switches change one property each:
#   --sigma full    L L' + diag(d): correlation plus genus-specific dispersion (the benchmark's setting)
#   --sigma factor  off-diagonal only (L L'): between-genus correlation, nothing else
#   --sigma 0       Sigma = 0: the plain Dirichlet-multinomial at the template's fitted A (no random draws for Z)
#   --vary_A 1      do not renormalise alpha_i to the fixed total A, so per-sample precision varies as well
#
# Template: <processed>/<template>/dirichlet_<control>.csv (alpha_dm_mle per genus, from estimate_dirichlet.R) and the
# control-group library sizes from <processed>/<template>/{Y.csv,X.csv}. All template genera are simulated.
#
# Strata by mean proportion over all genera: HIGH > 1%, MED 0.1-1%, LOW < 0.1%. Each case (simulation_cases.csv)
# draws n_high / n_med / n_low genera to perturb and applies, by direction:
#   up        alpha_j * exp(lfc_j), lfc_j ~ Uniform on the stratum's log fold-change range
#   balanced  the same magnitudes, half of each stratum up and half down (an odd one is randomised)
# Rare LOW genera perturbed upward are often absent from every control sample yet present in treated ones; these
# look like structural zeros to the methods and arise naturally, so no dedicated structural-zero case is simulated.
# Fold-change ranges by stratum: HIGH 1.25-2x, MED 1.5-4x, LOW 2-8x (log-uniform). Library sizes are resampled from
# the control samples. Only genera with zero counts in every sample of a replicate are dropped (as in the miRNA
# benchmark); no prevalence filter is applied -- methods are run as shipped, and a genus a method declines to test
# counts as not called (scripts/score_microbiome.R).
#
# Output (the layout nbsr-benchmarks reads): <out>/<case>/sample_<n>/rep<k>/{Y.csv, X.csv, fc.csv, settings.json}
# plus <out>/replicate_summary.csv. fc.csv keeps the pipeline's column names: "miRNA" is the feature key,
# alpha_null / alpha_alt the concentrations; true_effect = log(alpha_alt/alpha_null).
#
# Usage: Rscript data_processor/simulate_lndm.R --lognormal <lognormal_H_A38_q5.rds> [--sigma full] [--vary_A 0]
#          [--template crc_baxter_genus] [--control H] [--out DIR]
#          [--config data_processor/simulation_cases.csv] [--cases A,B,...] [--sample_sizes 40,80,160,320] [--reps 20] [--seed 1]
#          [--A <concentration>] [--depth_scale 1]   (the concentration sweep: --sigma 0 --A 38/100/400/... )
#          [--A_logsd 0]   sample-specific concentration A_i, lognormal around A with this log-sd (the per-sample precision experiment)
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
setwd(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
opt <- function(flag, default) { a <- commandArgs(TRUE); i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
template     <- opt("--template", "crc_baxter_genus")
control      <- opt("--control", "H")
ln_file      <- opt("--lognormal", NA)
sigma_mode   <- opt("--sigma", "full")
vary_A       <- as.integer(opt("--vary_A", "0")) == 1L
out_root     <- opt("--out", file.path("data/microbiome/simulated", paste0(template, "_lndm_", sigma_mode, if (vary_A) "_varyA" else "")))
config_file  <- opt("--config", "data_processor/simulation_cases.csv")
sample_sizes <- as.integer(strsplit(opt("--sample_sizes", "40,80,160,320"), ",")[[1]])
n_reps       <- as.integer(opt("--reps", "20"))
base_seed    <- as.integer(opt("--seed", "1"))
A_override   <- as.numeric(opt("--A", NA))           # Dirichlet concentration; default: the lognormal file's A (or the DM fit's with --sigma 0)
depth_scale  <- as.numeric(opt("--depth_scale", "1"))  # multiply the resampled library sizes
A_logsd      <- as.numeric(opt("--A_logsd", "0"))     # per-sample concentration: A_i = A exp(sigma e_i - sigma^2/2), e_i ~ N(0,1); 0 = one A for every sample

STRATA_CUTS <- c(HIGH = 0.01, MED = 0.001)                  # mean proportion: HIGH > 1%, MED 0.1-1%, LOW < 0.1%
FC_RANGE    <- list(HIGH = log(c(1.25, 2)), MED = log(c(1.5, 4)), LOW = log(c(2, 8)))
MIN_SAMPLES <- 1                                            # drop a genus only if it has zero counts in every sample

# ---- template ----------------------------------------------------------------------------------------
d   <- file.path("data/microbiome/processed", template)
est <- fread(file.path(d, paste0("dirichlet_", control, ".csv")))
Y   <- as.matrix(read.csv(file.path(d, "Y.csv"), row.names = 1, check.names = FALSE))
X   <- read.csv(file.path(d, "X.csv"), row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(Y), rownames(X)), all(est$feature %in% rownames(Y)))
lib_sizes <- colSums(Y[, X$trt == control])
alpha <- est$alpha_dm_mle; names(alpha) <- est$feature; A0 <- sum(alpha)
# ---- lognormal layer ---------------------------------------------------------------------------------
stopifnot(sigma_mode %in% c("0", "factor", "full"))
if (sigma_mode != "0") {
  if (is.na(ln_file)) stop("--lognormal <file from estimate_lognormal.R> is required unless --sigma 0")
  LN <- readRDS(ln_file); stopifnot(identical(LN$genera, names(alpha)))
  L_mat <- LN$L; d_vec <- if (sigma_mode == "full") LN$d else rep(0, length(alpha)); A_ln <- LN$A
  s2_vec <- rowSums(L_mat^2) + d_vec            # total log-scale variance per genus, for the mean-preserving centring
} else { L_mat <- NULL; d_vec <- rep(0, length(alpha)); A_ln <- A0 }
# A used for the Dirichlet layer: the lognormal file's A (a knob there), the DM total otherwise
A_use <- if (is.na(A_override)) A_ln else A_override; alpha <- alpha * (A_use / A0); A0 <- A_use
if (depth_scale != 1) lib_sizes <- round(lib_sizes * depth_scale)
draw_Z <- function(n) { K <- length(alpha); if (sigma_mode == "0") return(NULL)      # no RNG use: the plain Dirichlet-multinomial draws exactly
  Z <- matrix(rnorm(K * n, 0, rep(sqrt(d_vec), n)), K, n)
  if (!is.null(L_mat)) Z <- Z + L_mat %*% matrix(rnorm(ncol(L_mat) * n), ncol(L_mat), n); Z }
cat(sprintf("Lognormal layer: sigma=%s, vary_A=%s, A=%.1f%s\n", sigma_mode, vary_A, A0,
            if (sigma_mode != "0") sprintf(", rank %d, median diag %.2f", ncol(L_mat), median(d_vec)) else ""))
est[, stratum := fifelse(mean_prop > STRATA_CUTS["HIGH"], "HIGH", fifelse(mean_prop > STRATA_CUTS["MED"], "MED", "LOW"))]
strata_n <- est[, .N, by = stratum][order(-N)]
cat(sprintf("Template %s (%s): %d genera, A = %.2f, library size median %d [%d, %d]\n", template, control, nrow(est), A0,
            as.integer(median(lib_sizes)), min(lib_sizes), max(lib_sizes)))
cat("Strata (mean proportion over all genera):", paste(sprintf("%s=%d", strata_n$stratum, strata_n$N), collapse = ", "), "\n")

CASES <- fread(config_file)
# Each case's seed comes from its own seed_id, so running a subset of cases (--cases) or adding a row
# leaves every other case's replicates unchanged.
stopifnot("seed_id" %in% names(CASES), !anyDuplicated(CASES$seed_id))
cases <- strsplit(opt("--cases", paste(CASES$case, collapse = ",")), ",")[[1]]
CASES <- CASES[case %in% cases]
stopifnot(nrow(CASES) > 0, all(CASES$direction %in% c("up", "balanced")))
for (s in c("HIGH", "MED", "LOW")) {
  need <- CASES[[paste0("n_", tolower(s))]]; have <- est[stratum == s, .N]
  if (any(need > have)) stop(sprintf("Case asks for %d %s genera but the template has %d", max(need), s, have))
}

# ---- simulation --------------------------------------------------------------------------------------
rdirichlet_cols <- function(n, a) { g <- matrix(rgamma(length(a) * n, shape = a, rate = 1), length(a), n); sweep(g, 2, colSums(g), "/") }
normalize <- function(x) x / sum(x)

simulate_one <- function(cs, n, rep, seed) {
  set.seed(seed)
  # draw the perturbed genera, stratified
  de <- character(0); dir <- character(0)
  for (s in c("HIGH", "MED", "LOW")) {
    k <- cs[[paste0("n_", tolower(s))]]; if (k == 0) next
    picked <- sample(est[stratum == s, feature], k)
    sgn <- if (cs$direction == "up") rep("up", k) else sample(rep(c("up", "down"), length.out = k))   # odd one lands at random
    de <- c(de, picked); dir <- c(dir, sgn)
  }
  names(dir) <- de
  str_de <- est$stratum[match(de, est$feature)]
  lfc <- vapply(seq_along(de), function(i) runif(1, FC_RANGE[[str_de[i]]][1], FC_RANGE[[str_de[i]]][2]), numeric(1))
  lfc <- ifelse(dir == "down", -lfc, lfc)
  alpha_alt <- alpha; alpha_alt[de] <- alpha[de] * exp(lfc)
  # counts: n control ("null") then n treated ("alt") samples, library sizes resampled from the template controls
  ns <- sample(lib_sizes, 2 * n, replace = TRUE)
  # per-sample Dirichlet parameters: the group's alpha jittered by the lognormal layer (Z = 0 when sigma = 0)
  Zc <- draw_Z(n); Zt <- draw_Z(n)
  # centred so that E[exp(Z_j - s2_j/2)] = 1: each genus keeps its mean share of the Dirichlet parameter
  alpha_i <- function(a, Z) { if (is.null(Z)) return(matrix(a, length(a), n)); M <- exp(log(a) - s2_vec / 2 + Z); if (vary_A) M else sweep(M, 2, colSums(M), "/") * sum(a) }
  Ac <- alpha_i(alpha, Zc); At <- alpha_i(alpha_alt, Zt)
  if (A_logsd > 0) {   # sample-specific precision, mean A preserved: every sample's alpha_i rescaled to its own total
    scale_A <- function(M) { Ai <- exp(A_logsd * rnorm(ncol(M)) - A_logsd^2 / 2); sweep(M, 2, Ai, "*") }
    Ac <- scale_A(Ac); At <- scale_A(At)
  }
  P  <- cbind(sapply(seq_len(n), function(i) rdirichlet_cols(1, Ac[, i])), sapply(seq_len(n), function(i) rdirichlet_cols(1, At[, i])))
  Yc <- sapply(seq_len(2 * n), function(i) rmultinom(1, ns[i], P[, i]))
  rownames(Yc) <- names(alpha); colnames(Yc) <- sprintf("S%03d", seq_len(2 * n))
  trt <- rep(c("null", "alt"), each = n)
  retained <- rowSums(Yc > 0) >= MIN_SAMPLES
  # truth
  fc <- data.table(miRNA = names(alpha), alpha_null = alpha, alpha_alt = alpha_alt,
                   alpha_null_bar = normalize(alpha), alpha_alt_bar = normalize(alpha_alt))
  fc[, log_fc := log(alpha_alt_bar / alpha_null_bar)]          # change in proportion (relative scale)
  fc[, true_effect := log(alpha_alt / alpha_null)]            # change in absolute abundance
  fc[, is_de := miRNA %in% de]
  fc[, direction := fifelse(is_de, dir[miRNA], "none")]
  fc[, stratum := est$stratum[match(miRNA, est$feature)]]
  fc[, mean_prop_template := est$mean_prop[match(miRNA, est$feature)]]
  fc[, prevalence_template := est$prevalence[match(miRNA, est$feature)]]
  fc[, expected_reads_null := alpha_null_bar * median(lib_sizes)]   # at the template's median library size
  fc[, null_all_zero := rowSums(Yc[, trt == "null"]) == 0]   # absent from every control sample (looks like a structural zero)
  fc[, alt_all_zero  := rowSums(Yc[, trt == "alt"])  == 0]
  fc[, retained := retained[miRNA]]
  list(Y = Yc[retained, ], X = data.frame(trt = trt, row.names = colnames(Yc)), fc = fc,
       shift = log(sum(alpha_alt) / A0), n_retained = sum(retained), n_de_retained = sum(retained[de]),
       zero_frac = mean(Yc[retained, ] == 0), pivot = tail(names(alpha)[retained], 1),
       A_i_sd = sd(colSums(Ac)),        # within the control group: 0 unless --vary_A 1 or --A_logsd > 0
       n_de_null_all_zero = sum(fc[is_de == TRUE & retained == TRUE, null_all_zero]), lfc_mean = mean(abs(lfc)))
}

summary_rows <- list()
for (ci in seq_len(nrow(CASES))) {
  cs <- CASES[ci]
  for (n in sample_sizes) for (rep in seq_len(n_reps)) {
    seed <- base_seed * 1e7 + cs$seed_id * 1e5 + n * 1e2 + rep   # seed_id is fixed per case in the config
    sim <- simulate_one(cs, n, rep, seed)
    dir_out <- file.path(out_root, cs$case, sprintf("sample_%d", n), sprintf("rep%d", rep))
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    write.csv(sim$Y, file.path(dir_out, "Y.csv"))
    write.csv(sim$X, file.path(dir_out, "X.csv"))
    fwrite(sim$fc, file.path(dir_out, "fc.csv"))
    write_json(list(template = template, control = control, template_file = file.path(d, paste0("dirichlet_", control, ".csv")),
                    case = cs$case, direction = cs$direction, n_high = cs$n_high, n_med = cs$n_med, n_low = cs$n_low,
                    sample_size = n, replicate = rep, seed = seed,
                    strata_cuts_mean_prop = as.list(STRATA_CUTS), fc_range = lapply(FC_RANGE, exp), min_samples = MIN_SAMPLES,
                    A_null = A0, A_alt = A0 * exp(sim$shift), shift_log_total = sim$shift,
                    lognormal = list(sigma = sigma_mode, vary_A = vary_A, depth_scale = depth_scale, A_logsd = A_logsd, file = if (sigma_mode != "0") ln_file else NULL,
                                     rank = if (sigma_mode != "0") ncol(L_mat) else 0L, A_used = A0, A_i_sd = sim$A_i_sd),
                    n_features_template = length(alpha), n_features_retained = sim$n_retained, n_de_retained = sim$n_de_retained,
                    pivot_last_feature = sim$pivot),
               file.path(dir_out, "settings.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
    summary_rows[[length(summary_rows) + 1]] <- data.table(case = cs$case, sample_size = n, replicate = rep,
      n_de = cs$n_high + cs$n_med + cs$n_low, n_de_retained = sim$n_de_retained, shift = sim$shift,
      n_de_null_all_zero = sim$n_de_null_all_zero, n_retained = sim$n_retained, zero_frac = sim$zero_frac, A_i_sd = sim$A_i_sd,
      mean_abs_lfc = sim$lfc_mean, pivot = sim$pivot)
  }
}
summ <- rbindlist(summary_rows); fwrite(summ, file.path(out_root, "replicate_summary.csv"))
cat("\nPer case (over sample sizes and replicates):\n")
print(summ[, .(n_de = n_de[1], de_retained = round(mean(n_de_retained), 1), de_absent_in_controls = round(mean(n_de_null_all_zero), 1),
               shift_mean = round(mean(shift), 3),
               shift_range = sprintf("[%.2f, %.2f]", min(shift), max(shift)), features_retained = round(mean(n_retained)),
               zero_frac = round(mean(zero_frac), 3), mean_abs_lfc = round(mean(mean_abs_lfc), 2)), by = case])
cat("\nWritten to", out_root, ":", nrow(summ), "replicate directories\n")
