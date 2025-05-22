
# Good-Toulmin Estimator
good_toulmin <- function(freq_table, m) {
  # freq_table: named vector where names are frequencies and values are counts
  # m: # of future samples (=t*n, t, extrapolation factor (e.g., t = 1 for doubling the sample size)
  # n: # of past samples
  unseen <- rep(0, length(m))
  n = sum(freq_table)
  t = m / n
  for (f in names(freq_table)) {
    ff = as.numeric(f)
    unseen <- unseen - ((-t)^ff) * freq_table[f]
  }
  return(unseen) # Return real part in case of complex result
}

# Efron-Thisted Estimator
efron_thisted <- function(freq_table, m, max_terms = 10) {
  # freq_table: named vector where names are frequencies and values are counts
  # m: # of future samples (=t*n, t, extrapolation factor (e.g., t = 1 for doubling the sample size)
  # n: # of past samples
  unseen <- rep(0, length(m))
  n = sum(freq_table)
  t = m / n
  for (f in names(freq_table)) {
    ff = as.numeric(f)
    S = 1 - pbinom(ff, max_terms, 1/(1+t))
    tmp = ((-t)^(ff+1)) * S
    idx = t > max_terms
    tmp[idx] = 0
    unseen <- unseen + freq_table[f] * tmp
  }  
  return(unseen)
}

# Example usage
# Suppose we have observed frequencies: 118 species seen 1 time, 74 species seen 2 times, etc.
freq_table <- c(118,	74,	44,	24,	29,	22,	20,	19,	20,	15,	12,	14,	6, 1, 6)
# Corbet -FIsher butterfly data
names(freq_table) = as.character(1:length(freq_table))
print(sum(freq_table))

m = seq(1, sum(freq_table) * 1.3)
gt_estimate <- good_toulmin(freq_table, m)
et_estimate <- efron_thisted(freq_table, m)
et_estimate_20 <- efron_thisted(freq_table, m, 20)

plot(range(m), range(c(gt_estimate, et_estimate)), type='n', ylab='Estimated unseen species', xlab='# of future samples')
lines(m, gt_estimate)
lines(m, et_estimate, col='blue')
lines(m, et_estimate_20, col='green')
# cat("Good-Toulmin estimate:", gt_estimate, "\n")
# cat("Efron-Thisted estimate:", et_estimate, "\n")
