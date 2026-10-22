// nbsr with the dispersion fixed at values supplied as data: the oracle
// ablation. Every phi_ij is given (log_phi_fixed), so the only parameters
// are the mean model's beta. Under the Dirichlet-multinomial generator the
// true precision of feature j in group g is (A_g + 1) pi_j^(g) / (1 - pi_j^(g))
// with pi^(g) the group's composition and A_g its total concentration;
// build_stan_data() computes it from fc.csv. If the dispersion model's
// misspecification is what biases the fitted composition, this fit is
// unbiased; nothing else changes.
data {
  int<lower=1> N;
  int<lower=2> J;
  int<lower=1> P;
  matrix[N, P] X;
  vector<lower=0>[N] s;
  array[N, J] int<lower=0> Y;
  vector<lower=0>[P] sigma_beta2;
  matrix[N, J] log_phi_fixed;   // log precision of every count, given
}

parameters {
  matrix[P, J-1] z_beta;
}

transformed parameters {
  matrix[P, J-1] beta_raw;
  for (p in 1:P) beta_raw[p, ] = z_beta[p, ] * sqrt(sigma_beta2[p]);
}

model {
  matrix[N, J-1] eta_nonref = X * beta_raw;
  to_vector(z_beta) ~ std_normal();
  for (i in 1:N) {
    vector[J] eta_i;
    vector[J] log_pi_i;
    eta_i[1:(J-1)] = eta_nonref[i, ]';
    eta_i[J] = 0.0;
    log_pi_i = log_softmax(eta_i);
    for (j in 1:J)
      Y[i, j] ~ neg_binomial_2_log(log(s[i]) + log_pi_i[j], exp(log_phi_fixed[i, j]));
  }
}
