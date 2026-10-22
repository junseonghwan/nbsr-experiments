// nbsr with a switch on what the dispersion model's b_pi term is evaluated at.
// Identical to nbsr.stan when disp_pi = 1. An ablation of the dispersion
// structure; the paper's model stays in nbsr.stan.
data {
  int<lower=1> N;            // samples
  int<lower=2> J;            // features (miRNAs); category J is the reference
  int<lower=1> P;            // covariates for the mean model
  int<lower=0> Q;            // EXTERNAL covariates for the dispersion model

  matrix[N, P] X;            // mean-model design matrix
  matrix[N, Q] W;            // dispersion external covariates (no intercept)
  vector<lower=0>[N] s;      // size factors
  array[N, J] int<lower=0> Y;

  vector<lower=0>[P] sigma_beta2;
  vector<lower=0>[Q] sigma_b2;

  real b_pi_loc;
  real<lower=0> b_pi_scale;

  // What b_pi multiplies:
  //   0  nothing: log_phi_ij = b_0 + b_j + W_i . b_w          (b_pi keeps its prior)
  //   1  logit pi_ij, the sample's own fitted composition     (= nbsr.stan)
  //   2  logit pi0_j, the composition implied by the intercept row alone
  //      (the reference level's composition), so phi does not depend on
  //      the treatment coefficients
  //   3  logit_pi_obs, the observed pooled proportion passed as data, so phi
  //      does not depend on any mean-model parameter
  int<lower=0, upper=3> disp_pi;
  vector[J] logit_pi_obs;

  // Per-sample multiplier on the per-feature offset b_j: all ones gives the
  // usual b_j for every sample; the treatment indicator gives an offset that
  // applies to treatment samples only (no slack in the control group).
  vector<lower=0, upper=1>[N] bj_mult;

  // Optional per-feature treatment offset gamma_j * trt_i in the dispersion,
  // with its own hierarchy (sigma_gamma, zero-centred). use_gamma = 0 gives
  // zero-length parameters, so the model is unchanged when it is off.
  int<lower=0, upper=1> use_gamma;
  vector<lower=0, upper=1>[N] trt_ind;
}

parameters {
  matrix[P, J-1] z_beta;

  real b_0;
  real b_pi;
  real<lower=0> sigma_bj;
  vector[J] z_bj;
  vector[Q] z_b;
  array[use_gamma] real<lower=0> sigma_gamma;
  vector[use_gamma ? J : 0] z_gamma;
}

transformed parameters {
  matrix[P, J-1] beta_raw;
  for (p in 1:P) beta_raw[p,] = z_beta[p,] * sqrt(sigma_beta2[p]);

  vector[Q] b_w;
  for (q in 1:Q) b_w[q] = z_b[q] * sqrt(sigma_b2[q]);

  vector[J] b_j = sigma_bj * z_bj;
  vector[use_gamma ? J : 0] gamma;
  if (use_gamma) gamma = sigma_gamma[1] * z_gamma;
}

model {
  matrix[N, J-1] eta_nonref;
  vector[N] w_disp;
  vector[J] logit_pi0;       // used when disp_pi == 2

  to_vector(z_beta) ~ std_normal();
  to_vector(z_b) ~ std_normal();
  z_bj ~ std_normal();

  b_0      ~ normal(0, 1);
  b_pi     ~ normal(b_pi_loc, b_pi_scale);
  sigma_bj ~ normal(0, 0.5);
  if (use_gamma) {
    z_gamma ~ std_normal();
    sigma_gamma[1] ~ normal(0, 0.5);
  }

  eta_nonref = X * beta_raw;
  w_disp     = W * b_w;

  {
    vector[J] eta0;
    vector[J] log_pi0;
    eta0[1:(J-1)] = beta_raw[1, ]';   // intercept row: the reference level's composition
    eta0[J] = 0.0;
    log_pi0 = log_softmax(eta0);
    logit_pi0 = log_pi0 - log1m_exp(log_pi0);
  }

  for (i in 1:N) {
    vector[J] eta_i;
    vector[J] log_pi_i;
    vector[J] disp_term;

    eta_i[1:(J-1)] = eta_nonref[i, ]';
    eta_i[J] = 0.0;
    log_pi_i = log_softmax(eta_i);

    if (disp_pi == 1) {
      disp_term = b_pi * (log_pi_i - log1m_exp(log_pi_i));
    } else if (disp_pi == 2) {
      disp_term = b_pi * logit_pi0;
    } else if (disp_pi == 3) {
      disp_term = b_pi * logit_pi_obs;
    } else {
      disp_term = rep_vector(0.0, J);
    }

    for (j in 1:J) {
      real log_phi_ij = b_0 + bj_mult[i] * b_j[j] + disp_term[j] + w_disp[i]
                       + (use_gamma ? trt_ind[i] * gamma[j] : 0.0);
      Y[i, j] ~ neg_binomial_2_log(log(s[i]) + log_pi_i[j], exp(log_phi_ij));
    }
  }
}
