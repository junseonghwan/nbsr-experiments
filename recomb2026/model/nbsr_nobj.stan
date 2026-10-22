// nbsr.stan without the per-feature dispersion offsets b_j (and their scale
// sigma_bj): log_phi_ij = b_0 + b_pi * logit(pi_ij) + W_i . b_w. An ablation
// to measure what the hierarchical offsets contribute at n = 10, where their
// posterior sd equalled their prior scale. The data block is identical to
// nbsr.stan so the same builder and priors serve both.
data {
  int<lower=1> N;            // samples
  int<lower=2> J;            // features (miRNAs); category J is the reference
  int<lower=1> P;            // covariates for the mean model
  int<lower=0> Q;            // EXTERNAL covariates for the dispersion model
                             //   (intercept + log w_i terms; NOT log pi, which is internal)

  matrix[N, P] X;            // mean-model design matrix
  matrix[N, Q] W;            // dispersion external covariates (do not include intercept): [log w_{i1}, log w_{i2}, ...]
  vector<lower=0>[N] s;      // size factors
  array[N, J] int<lower=0> Y;

  vector<lower=0>[P] sigma_beta2;
  vector<lower=0>[Q] sigma_b2;

  // Prior on b_pi: normal(b_pi_loc, b_pi_scale). Every prior hyperparameter
  // is data so that a variant is a registry entry, not another model file;
  // the defaults live in PRIOR_DEFAULTS in scripts/common.R.
  real b_pi_loc;
  real<lower=0> b_pi_scale;
}

parameters {
  matrix[P, J-1] z_beta;   // mean-model coefficients (non-reference features)

  real b_0;
  real b_pi;
  vector[Q] z_b;           // coefficients on the external covariates in W
}

transformed parameters {
  matrix[P, J-1] beta_raw;
  for (p in 1:P) beta_raw[p,] = z_beta[p,] * sqrt(sigma_beta2[p]);
  
  vector[Q] b_w;
  for (q in 1:Q) b_w[q] = z_b[q] * sqrt(sigma_b2[q]);
}

model {
  matrix[N, J-1] eta_nonref;
  vector[N] w_disp;          // external part of the dispersion predictor: W_i . b_w

  // priors
  to_vector(z_beta) ~ std_normal();
  to_vector(z_b) ~ std_normal();

  b_0      ~ normal(0, 1);
  b_pi     ~ normal(b_pi_loc, b_pi_scale);

  // linear predictors
  eta_nonref = X * beta_raw;   // mean model
  w_disp     = W * b_w;        // external dispersion covariates (per sample), reused across j

  // likelihood
  for (i in 1:N) {
    vector[J] eta_i;
    vector[J] log_pi_i;
    vector[J] logit_pi_i;

    eta_i[1:(J-1)] = eta_nonref[i, ]';
    eta_i[J] = 0.0;                       // reference category fixed at 0

    log_pi_i = log_softmax(eta_i);        // internal, model-dependent composition
    logit_pi_i = log_pi_i - log1m_exp(log_pi_i);

    for (j in 1:J) {
      real log_phi_ij = b_0 + b_pi * logit_pi_i[j] + w_disp[i];

      Y[i, j] ~ neg_binomial_2_log(log(s[i]) + log_pi_i[j], exp(log_phi_ij));
    }
  }
}
