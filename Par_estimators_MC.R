# Poisson-Gamma model
# Usage:
#> freq_table<-spec_byhost_l[['Escherichia']] %>% table  %>%  table
#> freq_table<-spec_byhost_l[['Escherichia']] %>% unlist %>% table  %>%  table
#> par = PoissonGamma_MLE(freq_table)
#> FisherPoissonGamma(par[1], par[2], freq_table, 10)

logp <- function(k, alpha, beta) {
  binom <- lgamma(alpha + k) - (lgamma(k + 1) + lgamma(alpha))
  successes <- alpha * (log(beta) - log(beta + 1))
  failures <- k * (log(1) - log(beta + 1))
  return(binom + successes + failures)
}

#' Negative log-likelihood function for the beta-binomial model
#'
#' @param alpha Numeric
#' @param beta Numeric
#' @param f Integer vector. Frequencies of observed counts (f_k for k = 1, 2, ...).
#' @return Numeric. The negative log-likelihood value.
negloglikelihood <- function(alpha, beta, f) {
  if ((alpha < 0) | (beta < 0)){return(10^+308)} else {
    logL <- 0
    logp0 <- log1p(- (beta / (beta + 1)) ^ alpha)
    for (k in seq_along(f)) {
      logL <- logL + f[k] * (logp(k, alpha, beta) - logp0)
    }
    return(-logL)
  }
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

library(optimx)
PoissonGamma_MLE<-function(freq_table,
                           debug=FALSE){
  #optim finds minima, so we negate the fct.
  fn <- function(par){
    negloglikelihood(par[1], par[2], freq_table)
  }
  initial_guess<-as.numeric(sample.int(5, 2))
  res<-optimr(initial_guess, fn, method = "Nelder-Mead")
  if (debug){print(res)}
  ret = res$par
  names(ret) = c('alpha', 'beta')
  return(ret)
}

FisherPoissonGamma<-function(alpha, beta, freq_table, m){
  #From Efron-Thisted Paper, page 438, below eq. 3.3
  gamma = beta / (1 - beta)
  n = c(freq_table %*% as.numeric(names(freq_table)))
  t = m / n
  eta1 = freq_table[1] # number of species observed exactly one time
  if (alpha == 0){
    return (eta1 / gamma * log(1 + gamma * t))
  }else{
    return(-eta1 * ((1+ gamma * t) ^ (-alpha) -  1) / (gamma * alpha))
  }
}

FisherPoissonGammaWrapper<-function(freq_table, m){
  par<-PoissonGamma_MLE(freq_table)
  FisherPoissonGamma(par[1], par[2], freq_table, m)
}
