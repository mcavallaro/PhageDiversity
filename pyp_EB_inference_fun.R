#' Empirical Bayes approach for Pitman-Yor prior

#' Just for balocchi_likelihood and rising_factorial
#source("NP_estimators_MC.R")

#' rising factorial needs to be able to deal with u=0
#' for u=0, it is defined as 0 as it is the empty product
#' should we add an integer check?
rising_factorial2 <-function(a, u){
  ifelse(u>0,prod(a + 0:(u-1)),0)
  }

#' I think the balocchi likelihood looks slightly different
#' 
#' Important: Input to this is the vector M1,m,...,Mn,n (Mr,n)
#' which records the number of species with count r
#' so this can contain zeroes
#' 
balocchi_likelihood2 <-function(M, alpha, theta){
  # from page 21 of ?Balocchi's paper
  # alpha in 0,1
  # theta > -alpha
  n = length(M)
  Sum = sum(M)
  Factor1 = factorial(n) * rising_factorial2(theta / alpha, Sum) / rising_factorial2(theta, n)
# the second factor differs from MC's version  
  for (i in 1:n){
    Factor1 = Factor1 * (alpha * rising_factorial2(1 - alpha, i-1) / factorial(i))^M[i] / factorial(M[i])
  }
  return(Factor1)
}

#' Function to extract Mr,n from species counts
#' argument: species counts  

extractM <- function(counts){
n <- sum(counts) #No of sampled indiv
M <-rep(0,n) #at least count 0 in each class
#for each species, +1 the right class of M
for (i in seq(along=counts)) {M[counts[i]] <- M[counts[i]]+1} 
return(M)
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

uhat_pyp <- function(alpha,theta,m,M){
  no_species <- sum(M)
  n <- sum(seq(along=m)*M) #how many indiv's sampled?
  uhat <- no_species + theta/alpha
  uhat <- uhat * (rising_factorial2(alpha + theta + n,m)/rising_factorial2(theta + n,m) - 1)
 return(uhat)
} 
