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
  unseen = ifelse(unseen < 0, 0, unseen)
  return(unseen)
}

# Efron-Thisted Estimator
# this truncates and weighs
# the summands of the Good-Toulmin estimator
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
  # remove negatives:
  unseen = ifelse(unseen < 0, 0, unseen)
  return(unseen)
}
#' Orlitsky, Suresh and Wu (2016) modification 
#' of the Good-Toulmin estimator
#' this averages over a random cutoff of 
#' how many summands of GT are used
#' As discussed in OSW (Table 1), we will use
#' their preferred cutoff via random
#' distribution (smoothed binomial)
#' Bin(k,2/(2+m/n) with 
#' k=floor(0.5*log(x = (n*t^2/(t-1)),base = 3)) 
#' We will just use it for t>1 and refer to GT for t<=1
SGT <-  function(freq_table, m) {
  # freq_table: named vector where names are frequencies and values are counts
  # m: # of future samples (=t*n, t, extrapolation factor (e.g., t = 1 for doubling the sample size)
  # n: # of past samples
  unseen <- rep(0, length(m))
  n = c(freq_table %*% as.numeric(names(freq_table)))
  t = m / n
  for (l in seq(along=unseen)){
  smooth1 <- rep(1,length(freq_table)) #=1 for all entries for GT
  names(smooth1) <- names(freq_table)
  if (t[l] > 1){#binomial smoothing
  k_osw <- floor(0.5*log(x = (n*t[l]^2/(t[l]-1)),base = 3))
  for (f in names(freq_table)) {
    ff = as.integer(f)
    #below: observe P(L\geq i)=P(L > i-1) 
    smooth1[f] <- pbinom(ff-1,size = k_osw,
                         prob = 2/(t[l]+2),lower.tail = FALSE
                         )
      }
  }
  for (f in names(freq_table)) {
    ff = as.numeric(f)
    unseen[l] <- unseen[l] - ((-t[l])^ff) * smooth1[f] * freq_table[f]
  }}
  unseen = ifelse(unseen < 0, 0, unseen)
  return(unseen)
}
