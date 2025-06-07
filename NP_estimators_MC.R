
# Good-Toulmin Estimator
good_toulmin <- function(freq_table, m) {
  # freq_table: named vector where names are frequencies and values are counts
  # m: # of future samples (=t*n, t, extrapolation factor (e.g., t = 1 for doubling the sample size)
  # n: # of past samples
  unseen <- rep(0, length(m))
  n = c(freq_table %*% as.numeric(names(freq_table)))
  t = m / n
  for (f in names(freq_table)) {
    ff = as.numeric(f)
    unseen <- unseen - ((-t)^ff) * freq_table[f]
  }
  return(unseen)
}

# Efron-Thisted Estimator
efron_thisted <- function(freq_table, m, max_terms = 10) {
  # freq_table: named vector where names are frequencies and values are counts
  # m: # of future samples (=t*n, t, extrapolation factor (e.g., t = 1 for doubling the sample size)
  # n: # of past samples
  unseen <- rep(0, length(m))
  n = c(freq_table %*% as.numeric(names(freq_table)))
  t = m / n
  for (f in names(freq_table)) {
    ff = as.numeric(f)
    S = 1 - pbinom(ff, max_terms, 1/(1+t))
    #tmp = ((-t)^(ff+1)) * S
    tmp = ((-t)^(ff)) * S
    idx = t > max_terms
    tmp[idx] = 0
    # unseen <- unseen + freq_table[f] * tmp
    unseen <- unseen - freq_table[f] * tmp
  }  
  return(unseen)
}

rising_factorial<-function(a, u){
  prod(a + 0:(u-1))
}

balocchi_likelihood<-function(M, alpha, theta){
  # from page 21 of ?Balocchi's paper
  # alpha in 0,1
  # theta > -alpha
  n = length(M)
  Sum = sum(M)
  Factor1 = factorial(n) * rising_factorial(theta / alpha, Sum) / rising_factorial(theta, n)
  for (i in 1:n){
    Factor1 = Factor1 * (rising_factorial(alpha * (1 - alpha), i-1) / factorial(i))^M[i] / factorial(M[i])
  }
  return(Factor1)
}


#Poisson-Gamma model

logp <- function(k, params) {
  alpha <- params[1]
  beta <- params[2]
  binom <- lgamma(alpha + k) - (lgamma(k + 1) + lgamma(alpha))
  successes <- alpha * (log(beta) - log(beta + 1))
  failures <- k * (log(1) - log(beta + 1))
  return(binom + successes + failures)
}

#' Negative log-likelihood function for the beta-binomial model
#'
#' @param params Numeric vector of length 2. Parameters alpha and beta.
#' @param f Integer vector. Frequencies of observed counts (f_k for k = 1, 2, ...).
#' @return Numeric. The negative log-likelihood value.
negloglikelihood <- function(params, f) {
  logL <- 0
  logp0 <- log1p(- (params[2] / (params[2] + 1)) ^ params[1])
  for (k in seq_along(f)) {
    logL <- logL + f[k] * (logp(k, params) - logp0)
  }
  return(-logL)
}

#' Estimate the number unobserved species using the beta-binomial model
#'
#' @param f Integer vector. Frequencies of observed counts (f_k for k = 1, 2, ...). 
#' frequency counts (value in row k corresponds to the number of clones with exactly k individuals in the sample)
#' @return Numeric vector. Estimated parameters alpha and beta, and total number of unobserved species.
# estimate <- function(f) {
#   fn <- function(par){negloglikelihood(par, f)}
#   res <- optim(par = c(50, 1), fn = fn)
#   alpha <- res$par[1]
#   beta <- res$par[2]
#   p0 <- (beta / (beta + 1)) ^ alpha
#   Sobs <- sum(f)
#   return(alpha, beta, Sobs / (1 - p0) - Sobs)
# }


FisherPoissonGamma<-function(freq_table, m, alpha, beta){
  #From BredleyEfron Paper, page 438
  
  fn <- function(par){negloglikelihood(par, f)}
  res <- optim(par = c(50, 1), fn = fn)
  alpha <- res$par[1]
  beta <- res$par[2]
  
  gamma = beta / (1 - beta)
  n = c(freq_table %*% as.numeric(names(freq_table)))
  t = m / n
  eta1 = f[1] # number ofspeciesobservedexactly 1 time
  if (alpha == 0){
    return (eta1 / gamma * log(1 + gamma * t))
  }else{
    return(-eta1 * ((1+ gamma * t) ^ (-alpha) -  1) (gamma * alpha))
  }
}


