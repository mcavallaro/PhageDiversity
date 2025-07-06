simulate_P <- function(alpha, theta, n = 1000) {
  #n is number of species
  if (alpha < 0 || alpha >= 1) stop("alpha must be in [0,1)")
  if (theta <= -alpha) stop("theta must be > -alpha")
  
  v <- numeric(n)
  w <- numeric(n)
  prod_term <- 1
  
  for (i in 1:n) {
    v[i] <- rbeta(1, 1 - alpha, theta + i * alpha)
    w[i] <- v[i] * prod_term
    prod_term <- prod_term * (1 - v[i])
  }
  
  return(w)
}


simulate_pitman_yor <- function(alpha, theta, base_distribution, n = 1000) {
  weights <- simulate_P(alpha, theta, n)
  S <- base_distribution(n)  # Sample from base distribution (e.g., atomic or non atomic??) (eta)

  Samples<-sample(n,n,replace = T, prob=weights)
  
  return(list(S[Samples], Samples))
  
}

# Example usage
alpha <- 0.5
theta <- 1.0
weights <- simulate_P(alpha, theta, 100)

barplot(weights[1:30], main="weights", ylab="Weight")


py_result <- simulate_pitman_yor(alpha = 0.5, theta = 1.0, base_distribution=rnorm, n = 100)
plot(table(py_result[[1]]))



source("utils.R")
simulatePoissonGamma<-function(shape, rate, n = 1000){
  
  lambdas<-rgamma(shape, rate, n = n)
  # shape is alpha
  # rate is beta
  # n is number of species
  
  # Draw Poisson samples with these lambdas
  speccounts <-rpois(lambda = lambdas, n = n)
  cat("The simulated total number of individuals is: ", sum(speccounts), "\n")
  freq_table<-getFrequencyTable(speccounts)
  return(freq_table)
}

