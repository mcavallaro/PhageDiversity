#Example usage
# Suppose we have observed frequencies: 118 species seen 1 time, 74 species seen 2 times, etc.
freq_table <- c(118,	74,	44,	24,	29,	22,	20,	19,	20,	15,	12,	14,	6, 1, 6)
# Corbet -FIsher butterfly data
names(freq_table) = as.character(1:length(freq_table))


cat("\n", n1)
cat(" n of species: ", sum(freq_table))
cat(". n of isolates: ", c(freq_table %*% as.numeric(names(freq_table))))

n_isolates =  c(freq_table %*% as.numeric(names(freq_table)))

source("NP_estimators_MC.R")


m = seq(1, n_isolates * 1.3)
gt_estimate <- good_toulmin(freq_table, m)
et_estimate <- efron_thisted(freq_table, m)
et_estimate_20 <- efron_thisted(freq_table, m, 20)

new_species_SGT = vector(mode = "numeric", length = length(m))
for (i in m){
  new_species_SGT[i]=sgt_Delta(r = as.integer(names(freq_table)),
                               N_r = freq_table,
                               m=n_isolates, t=i/n_isolates, adj = F)
}


plot(range(m), range(c( et_estimate)), type='n', main='Butterfly data',
     ylab='Estimated unseen species', xlab='# of future samples')

points(m, new_species_SGT)
abline(v=n_isolates, col='red')
lines(m, gt_estimate)
lines(m, et_estimate, col='blue')
lines(m, et_estimate_20, col='green')


legend('topleft', lty=c(NA,1,1,1,1), pch=c(1,NA,NA,NA,NA), c('SGT_delta','GT', 'ET', 'ET (max_terms=20)', 't=1'),
       col=c('black','black', 'blue', 'green', 'red'))


# cat("Good-Toulmin estimate:", gt_estimate, "\n")
# cat("Efron-Thisted estimate:", et_estimate, "\n")
