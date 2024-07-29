digamma_inv <- function(y, gam = -digamma(1)) {
  x <- exp(y) + 0.5
  idx <- (y < -2.22)
  x[idx] <- -1/(y[idx] + gam)
  return(x)
}

log_dirichlet <- function(p, alpha)
{
  alpha0 <- sum(alpha)
  log_lik <- lgamma(alpha0) - sum(lgamma(alpha)) + sum((alpha - 1) * log(p))
  return(log_lik)
}
optimize_concentration <- function(Y, alpha_bar, s_init, max_iter = 1000)
{
  s_old <- s_init
  nn <- colSums(Y)
  for (i in 1:max_iter)
  {
    ret <- alpha_bar * digamma(Y + s_old * alpha_bar)
    num <- sum(alpha_bar * digamma(Y + s_old * alpha_bar) - alpha_bar * digamma(s_old * alpha_bar))
    den <- sum(digamma(nn + s_old) - digamma(s_old))
    s_new <- s_old * num / den
    #print(paste0(s_old, "->", s_new))
    if (abs(s_new - s_old) < 1e-6) {
      print("Converged.")
      break
    }
    s_old <- s_new
  }
  return(list(s=s_old, iters=i))
}

newton <- function(alpha_old, y, iter=10)
{
  for (i in 1:iter)
  {
    alpha_new <- alpha_old - (digamma(alpha_old) - y) / trigamma(alpha_old)
    alpha_old <- alpha_new  
    #print(sum(alpha_new))
  }
  return(alpha_new)
}

estimate_dirichlet <- function(Y, initial_alpha, max_iters = 50, verbose = FALSE)
{
  # Smooth it by adding a 1.
  Y <- Y + 1e-6
  props <- apply(Y, 1, function(x) {
    x / sum(x)
  })

  # Compute geometric mean of each category.
  log_geom_mean <- rowMeans(log(props))
  
  alpha_old <- initial_alpha
  y <- digamma(sum(alpha_old)) + log_geom_mean
  alpha_old <- digamma_inv(y)
  s_prev <- sum(alpha_old)
  converged <- FALSE

  for (j in 1:max_iters)
  {
    y <- digamma(sum(alpha_old)) + log_geom_mean
    alpha_new <- newton(alpha_old, y, 20)
    alpha_old <- alpha_new
    s <- sum(alpha_new)
    if (verbose) {
      print(s)
    }
    if (abs(s - s_prev) < 1e-10) {
      converged <- TRUE
      print("Converged")
      break
    }
    s_prev <- s
  }
  
  alpha_est1 <- alpha_new
  #m_est1 <- alpha_est1 / sum(alpha_est1)
  return(list(alpha_hat=alpha_est1, converged = converged))
}
