# The estimand, shared by every inference method.
#
# Delta_adj is the per-miRNA treatment effect in CLR space, with each draw
# re-centred on the mode of that draw's effect distribution, which is what
# makes effects comparable across draws under a compositional constraint.
#
# This is a property of the model, not of the algorithm that fit it. Keeping it
# in one place is what guarantees NBSR-HMC, NBSR-Laplace and any future method
# are scored by identical code rather than by copies that drift apart.

# Column of the design matrix holding the treatment indicator (prepare_inputs
# builds it as [intercept, trt]).
TREATMENT_COL <- 2

# Mode of a Gaussian KDE with the default bandwidth. density() evaluates a
# binned approximation on a 512-point grid, so its argmax is only the nearest
# grid point; maximizing the exact KDE within the neighbouring cells returns
# the mode itself, removing a quantization error of up to one grid step
# (~range/500).
kde_mode <- function(v) {
  kde <- density(v)
  k <- which.max(kde$y)
  lo <- kde$x[max(k - 1L, 1L)]
  hi <- kde$x[min(k + 1L, length(kde$x))]
  optimize(function(m) sum(dnorm(m, v, kde$bw)), c(lo, hi), maximum = TRUE)$maximum
}

# The centring locates the null miRNAs' common compositional shift as the
# mode of the effect distribution over a reference set of miRNAs. The set is
# the most abundant fraction of the panel, ranked by pooled proportion so
# that neither condition nor truth is involved. Those effects are estimated
# precisely and, under the model's own dispersion misspecification, without
# bias; null miRNAs of intermediate abundance carry a small positive bias
# that grows sharper with n, and a mode taken over the whole panel drifts
# toward them and passes the drift to the perturbed miRNAs as an offset of
# about -0.03 (0.87 coverage on HIGH at n = 40). The mode is robust to the
# perturbed miRNAs that land in a small set; a mean is not. Settable through
# NBSR_CENTRING_FRACTION; 1 restores the whole panel. The paper's section
# "Where the shift is estimated from" gives the sensitivity table.
CENTRING_FRACTION_DEFAULT <- 0.1
CENTRING_MIN_MIRNAS <- 30

centring_fraction <- function() {
  from_env <- Sys.getenv("NBSR_CENTRING_FRACTION", unset = "")
  if (!nzchar(from_env)) return(CENTRING_FRACTION_DEFAULT)
  value <- suppressWarnings(as.numeric(from_env))
  if (is.na(value) || value <= 0 || value > 1) {
    stop("NBSR_CENTRING_FRACTION must be in (0, 1]; got: ", from_env, call. = FALSE)
  }
  value
}

# Indices of the reference set: the top `fraction` of miRNAs by abundance,
# never fewer than CENTRING_MIN_MIRNAS (or all of them, if the panel is
# smaller than that).
centring_set <- function(abundance, fraction) {
  J <- length(abundance)
  k <- max(min(J, CENTRING_MIN_MIRNAS), ceiling(fraction * J))
  order(abundance, decreasing = TRUE)[seq_len(k)]
}

# beta_draws is (draws, P, J - 1) in ALR space, omitting the reference
# category. `abundance` ranks the J miRNAs (in the same order as the model's
# categories, reference last) for the centring set. Returns Delta_adj
# (draws, J), the mode-centred treatment effect. Only that column is kept:
# the full (draws, J, P) CLR array is P times larger and nothing downstream
# reads it.
delta_from_beta <- function(beta_draws, abundance = NULL,
                            fraction = centring_fraction()) {
  S <- dim(beta_draws)[1]
  P <- dim(beta_draws)[2]
  J <- dim(beta_draws)[3] + 1   # reference category is added back below
  if (is.null(abundance)) {
    if (fraction < 1) stop("A centring fraction below 1 needs the miRNA abundances.", call. = FALSE)
    ref <- seq_len(J)
  } else {
    if (length(abundance) != J) stop("abundance must have one value per miRNA (", J, ").", call. = FALSE)
    ref <- centring_set(abundance, fraction)
  }

  # The estimand is column TREATMENT_COL of the design. With P == 1 there is
  # no such column, and [i, , ] would also collapse to a vector, so the loop
  # below would either index out of bounds or build the wrong shape.
  if (P < 2) {
    stop("beta_draws has ", P, " covariate(s); the treatment effect is column ",
         TREATMENT_COL, " of the design.", call. = FALSE)
  }

  Delta_adj <- array(NA_real_, dim = c(S, J))

  for (i in seq_len(S)) {
    b <- rbind(t(beta_draws[i, , ]), 0)          # add the reference category
    x <- scale(b, center = TRUE, scale = FALSE)  # ALR -> CLR
    Delta_adj[i, ] <- x[, TREATMENT_COL] - kde_mode(x[ref, TREATMENT_COL])
  }

  list(Delta_adj = Delta_adj, n_draws = S, centring_fraction = fraction,
       centring_set_size = length(ref))
}

# Convergence diagnostics on the estimand itself rather than on the raw
# parameters, since Delta_adj is what the results actually rest on.
delta_adj_diagnostics <- function(Delta_adj, num_chains) {
  S <- nrow(Delta_adj)
  if (S %% num_chains != 0) {
    stop("Draw count ", S, " is not divisible by ", num_chains,
         " chains; cannot reconstruct the chains for diagnostics.", call. = FALSE)
  }
  num_iter <- S %/% num_chains
  J <- ncol(Delta_adj)

  arr <- array(
    Delta_adj,
    dim = c(num_iter, num_chains, J),
    dimnames = list(
      iteration = seq_len(num_iter),
      chain     = paste0("chain", seq_len(num_chains)),
      variable  = sprintf("Delta_adj[%d]", seq_len(J))
    )
  )
  # Named so the columns come out as rhat/ess_bulk/... rather than being
  # labelled with the namespaced expression that produced them.
  posterior::summarise_draws(
    posterior::as_draws_array(arr),
    rhat      = posterior::rhat,
    ess_bulk  = posterior::ess_bulk,
    ess_tail  = posterior::ess_tail,
    mcse_mean = posterior::mcse_mean
  )
}

SUMMARY_PROBS <- c(q005 = 0.005, q025 = 0.025, q05 = 0.05, q10 = 0.10,
                   q25 = 0.25, q50 = 0.50, q75 = 0.75, q90 = 0.90,
                   q95 = 0.95, q975 = 0.975, q995 = 0.995)

# Reduce Delta_adj to one row per miRNA, tagged with the method and the cell of
# the simulation grid it came from, then joined to the true effects. Carrying
# method_id in the table itself is what makes a summary self-describing, so
# combining across methods needs no path parsing.
summarize_delta <- function(Delta_adj, mirna, method_id, simul_case, sample_size,
                            replicate_id, true_fc, centring_fraction = 1) {
  qmat <- matrixStats::colQuantiles(Delta_adj, probs = SUMMARY_PROBS)
  colnames(qmat) <- names(SUMMARY_PROBS)

  dt <- data.table::data.table(
    miRNA        = mirna,
    method_id    = method_id,
    simul_case   = simul_case,
    sample_size  = sample_size,
    replicate_id = replicate_id,
    post_mean    = matrixStats::colMeans2(Delta_adj),
    post_sd      = matrixStats::colSds(Delta_adj),
    p_pos        = colMeans(Delta_adj > 0)
  )
  dt <- cbind(dt, qmat)
  dt[, centring_fraction := centring_fraction]
  dt <- merge(dt, true_fc, by = "miRNA", all.x = TRUE, all.y = FALSE)

  # The target the estimand is scored against. fc.csv's log_fc is the change
  # in *proportion*, log(alpha_alt_bar / alpha_null_bar): when a few miRNAs go
  # up, every other miRNA's proportion falls by a common factor, so log_fc is
  # non-zero for miRNAs whose abundance never changed. Mode-centring Delta_adj
  # removes exactly that common factor, so the estimand is on the *absolute*
  # scale, log(alpha_alt / alpha_null): zero for every non-DE miRNA. Carrying
  # that column here means no downstream analysis can score against log_fc
  # by mistake.
  dt[, true_effect := log(alpha_alt / alpha_null)]
  dt[, is_de := abs(true_effect) > 1e-12]
  dt
}
