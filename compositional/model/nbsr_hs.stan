// nbsr with a regularised horseshoe prior on the treatment row of beta.
//
// Identical to nbsr.stan except for the prior on beta_raw[2, ]: instead of
// normal(0, sqrt(sigma_beta2[2])) the treatment effects get
//   beta_raw[2, j] = tau * lambda_tilde_j * z_beta[2, j]
// with local scales lambda_j ~ half-Cauchy(0, 1), global scale
// tau ~ half-Cauchy(0, hs_tau0), and the slab regularisation of Piironen
// and Vehtari (2017), lambda_tilde^2 = c^2 lambda^2 / (c^2 + tau^2 lambda^2),
// c^2 = hs_slab_scale^2 * caux, caux ~ inv-gamma(nu/2, nu/2). The intercept
// row (and any further rows) keep the normal prior. hs_tau0 encodes the
// prior guess of the number of non-null effects: p0 / (K - p0) times a
// noise scale.
data {
  int<lower=1> N;
  int<lower=2> J;
  int<lower=2> P;            // row 2 is the treatment indicator
  int<lower=0> Q;

  matrix[N, P] X;
  matrix[N, Q] W;
  vector<lower=0>[N] s;
  array[N, J] int<lower=0> Y;

  vector<lower=0>[P] sigma_beta2;   // used for every row except row 2
  vector<lower=0>[Q] sigma_b2;
  real b_pi_loc;
  real<lower=0> b_pi_scale;

  real<lower=0> hs_tau0;
  real<lower=0> hs_slab_scale;
  real<lower=0> hs_slab_df;
}

parameters {
  matrix[P, J-1] z_beta;
  vector<lower=0>[J-1] lambda;
  real<lower=0> tau;
  real<lower=0> caux;

  real b_0;
  real b_pi;
  real<lower=0> sigma_bj;
  vector[J] z_bj;
  vector[Q] z_b;
}

transformed parameters {
  real<lower=0> c = hs_slab_scale * sqrt(caux);
  vector<lower=0>[J-1] lambda_tilde = sqrt(square(c) * square(lambda) ./ (square(c) + square(tau) * square(lambda)));
  matrix[P, J-1] beta_raw;
  for (p in 1:P) {
    if (p == 2) beta_raw[p, ] = tau * (z_beta[p, ] .* lambda_tilde');
    else beta_raw[p, ] = z_beta[p, ] * sqrt(sigma_beta2[p]);
  }

  vector[Q] b_w;
  for (q in 1:Q) b_w[q] = z_b[q] * sqrt(sigma_b2[q]);

  vector[J] b_j = sigma_bj * z_bj;
}

model {
  matrix[N, J-1] eta_nonref;
  vector[N] w_disp;

  to_vector(z_beta) ~ std_normal();
  lambda ~ cauchy(0, 1);
  tau    ~ cauchy(0, hs_tau0);
  caux   ~ inv_gamma(0.5 * hs_slab_df, 0.5 * hs_slab_df);
  to_vector(z_b) ~ std_normal();
  z_bj ~ std_normal();

  b_0      ~ normal(0, 1);
  b_pi     ~ normal(b_pi_loc, b_pi_scale);
  sigma_bj ~ normal(0, 0.5);

  eta_nonref = X * beta_raw;
  w_disp     = W * b_w;

  for (i in 1:N) {
    vector[J] eta_i;
    vector[J] log_pi_i;
    vector[J] logit_pi_i;

    eta_i[1:(J-1)] = eta_nonref[i, ]';
    eta_i[J] = 0.0;
    log_pi_i = log_softmax(eta_i);
    logit_pi_i = log_pi_i - log1m_exp(log_pi_i);

    for (j in 1:J) {
      real log_phi_ij = b_0 + b_j[j] + b_pi * logit_pi_i[j] + w_disp[i];
      Y[i, j] ~ neg_binomial_2_log(log(s[i]) + log_pi_i[j], exp(log_phi_ij));
    }
  }
}
