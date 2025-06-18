#' Empirical Bayes approach for Pitman-Yor prior
# Usage:
#> M = spec_byhost_l[['Escherichia']] %>% unlist() %>% table  %>% extractM()
#> par = PYP_MLE(M)
#> uhat_pyp(par[1], par[2], M, 10)



#' Just for balocchi_likelihood and rising_factorial
#source("NP_estimators_MC.R")

#' rising factorial needs to be able to deal with u=0
#' for u=0, it is defined as 0 as it is the empty product
#' should we add an integer check?
rising_factorial <-function(a, u){
  ifelse(u>0, prod(a + 0:(u-1)), 1)
}

log_rising_factorial <-function(a, u){ # I think we needed the log inside the sum ;D
  ifelse(u>0, sum(log(a + 0:(u-1))), 0)
}

#' I think the balocchi likelihood looks slightly different
#' 
#' Important: Input to this is the vector M1,m,...,Mn,n (Mr,n)
#' which records the number of species with count r
#' so this can contain zeroes
#' 
balocchi_likelihood <-function(M, alpha, theta, N){
  # from page 21 of ?Balocchi's paper
  # alpha in 0,1
  # theta > -alpha
  n = length(M)
  Sum = sum(M)
  # Factor1 = factorial(n) * rising_factorial2(theta / alpha, Sum) / rising_factorial2(theta, n)
  Factor1 = rising_factorial2(theta / alpha, Sum) / rising_factorial2(theta, n) # Factorial not needed for MLE
# the second factor differs from MC's version  
  cat('rrrr', rising_factorial2(theta / alpha, Sum), "\n")
  cat('rrrr', rising_factorial2(theta, n), "\n")  
  for (i in 1:N){
    Factor1 = Factor1 * (alpha * rising_factorial2(1 - alpha, i-1) / factorial(i))^M[i] / factorial(M[i])
    cat(i, Factor1, "\n")
  }
  return(Factor1)
}


neglog_balocchi_likelihood <-function(M, alpha, theta, N=NULL){
  if (!(alpha <=1 & alpha >=0 & theta > -alpha)){return(10^+308)} else {
  # from page 21 of ?Balocchi's paper
  # alpha in 0,1
  # theta > -alpha
  # usage log_balocchi_likelihood(M, 0.5, 1, 20) %>% exp()
  n = length(M)
  Sum = sum(M)
  if(is.null(N)){
    N = n
  }
  # Factor1 = factorial(n) * rising_factorial2(theta / alpha, Sum) / rising_factorial2(theta, n)
  Term = log_rising_factorial(theta / alpha, Sum) - log_rising_factorial(theta, n) # Factorial not needed for MLE
  # the second factor differs from MC's version  
  for (i in 1:N){
    # Term1 = Term1 +  (alpha * rising_factorial2(1 - alpha, i-1) / factorial(i))^M[i] / factorial(M[i])
    Term = Term + ( M[i] * (log(alpha) + log_rising_factorial(1 - alpha, i-1) - lfactorial(i) ) - lfactorial(M[i]) )
  }}
  return(-Term)
}


#' Function to extract Mr,n from species counts
#' argument: species counts  
extractM <- function(counts){
  #M = spec_byhost_l[['Escherichia']] %>% unlist() %>% table  %>% extractM()
  n <- sum(counts) #No of sampled indiv
  M <-rep(0,n) #at least count 0 in each class
  #for each species, +1 the right class of M
  for (i in seq(along=counts)) {M[counts[i]] <- M[counts[i]]+1} 
  return(M)
}

library(optimx)
PYP_MLE<-function(M){
  #optim finds minima, so we negate the fct.
  fn <- function(par){
    neglog_balocchi_likelihood(M, par[1], par[2])
  }
  # res <- optim(par = c(0.5, 1), fn = fn)
  res<-optimr(c(0.5, 1), fn, method = "Nelder-Mead") #, lower = c(0, -1), upper=c(1, 1000))
  print(res)
  ret = res$par
  names(ret) = c('alpha', 'theta')
  return(ret)
}

#' We will use the Empirical Bayes procedure
#' Step 1: ML estimate of PYP parameters, optimisation as in 
#' https://github.com/cecilia-balocchi/OrderedSSP
#'
#' TO DO
 
#' Step 2: Use Eq. 19, 20 from Balocchi paper 1 to estimate 
#' unseen species properties from the empirical parameter estimate
#' for now, estimate mean number of new species expected under
#' best estimate of alpha, theta
#' M is the vector M1,n,...Mn,n of species with specific count 1,...,n
uhat_pyp <- function(alpha, theta, M, m){
  no_species <- sum(M)
  n <- sum(seq(along=m)*M) #how many indiv's sampled?
  uhat <- no_species + theta/alpha
  uhat <- uhat * (rising_factorial(alpha + theta + n,m)/rising_factorial(theta + n,m) - 1)
 return(uhat)
}
