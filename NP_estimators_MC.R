
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

