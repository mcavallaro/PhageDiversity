#Example usage
# Suppose we have observed frequencies: 118 species seen 1 time, 74 species seen 2 times, etc.
freq_table <- c(118,	74,	44,	24,	29,	22,	20,	19,	20,	15,	12,	14,	6, 1, 6)
# Corbet -FIsher butterfly data
names(freq_table) = as.character(1:length(freq_table))
print(sum(freq_table))

source("NP_estimators_MC.R")


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
