library(magrittr)
library(dplyr)

# import functions for non-paramteric estimates
source("NP_estimators_MC.R")

# import functions for paramteric estimates
source("Par_estimators_MC.R")

# import functions for PYP estimates
source("pyp_EB_inference_fun.R")

# import functions for bootstrap and other utils
source("utils.R")

# import data
fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost <- fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host
#' count number of samples per host species
#' which ones are > 1000?
nosamples_host <- sapply(spec_byhost_l, nrow)
nospecies_host <- sapply(spec_byhost_l, function(x){nrow(unique(x))})
samples_1K <- which(nosamples_host>999)

# At this point, for each host, we have a collection of phages (individual observations, along with their species names)
# > spec_byhost_l[["Pseudomonas"]]
# # A tibble: 1,644 × 1
# `Phage Species`
# <chr>          
#   1 vOTU_0044      
# 2 vOTU_0044
# 3 vOTU_0045
# 4 vOTU_0048
# 5 vOTU_0049
# 6 vOTU_0048
# 7 vOTU_0049
# 8 vOTU_0052
# 9 vOTU_0054
# 10 vOTU_0056
# # ℹ 1,634 more rows
# # ℹ Use `print(n = ...)` to see more rows


# Estimate the missing species from the Sept2024 data and see how it compares with 
# the last dataset
new_species_GT = list()
new_species_ET = list()
new_species_PoiGamma = list()
new_species_PYP = list()

num_boostrap_samples = 50

for (n1 in names(samples_1K[c(1,2,3,4)])){
  # speccounts <- as.vector(unname(table(spec_byhost_l[[n1]])))
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  cat("\n", n1)
  cat(". n of species: ", sum(freq_table), " ", length(speccounts))
  m =  c(freq_table %*% as.numeric(names(freq_table)))
  cat(". n of isolates: ", m, " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
  m = seq(1, m * 1.3)
  new_species_GT[[n1]] = good_toulmin(freq_table, m)
  new_species_ET[[n1]] = efron_thisted(freq_table, m)
  new_species_PoiGamma[[n1]] = FisherPoissonGammaWrapper(freq_table, m)
  # M<-extractM(speccounts)
  # new_species_PYP[[n1]] = BalocchiPYPWrapper(M, m)
  
  tmp=bootstrapSpeciesCount(speccounts, num_boostrap_samples)
  for(i in 1:num_boostrap_samples){
    sample = tmp[[i]]
    freq_table<-getFrequencyTable(sample)
    # cat("\n", n1, " Bootstrap sample ", i)
    # cat(". n of species: ", sum(freq_table), " ", length(speccounts))
    # cat(". n of isolates: ", c(freq_table %*% as.numeric(names(freq_table))), " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
    m =  c(freq_table %*% as.numeric(names(freq_table)))
    m = seq(1, m * 1.3)
    attribute_name = sprintf("Boostrap sample %d", i)
    attr(new_species_GT[[n1]], attribute_name) = good_toulmin(freq_table, m)
    attr(new_species_ET[[n1]], attribute_name) = efron_thisted(freq_table, m)
    attr(new_species_PoiGamma[[n1]], attribute_name) = FisherPoissonGammaWrapper(freq_table, m)
  }
}


# plot
tab.blue="#1f77b4"
tab.orange="#ff7f0e"
tab.green="#2ca02c"
tab.red="#d62728"
tab.purple="#9467bd"
tab.blue.alpha="#1f77b444"
tab.orange.alpha="#ff7f0e44"
tab.green.alpha="#2ca02c44"
tab.red.alpha="#d6272844"
tab.purple.alpha="#9467bd44"
BootstrapPredictionIntervals<-function(new_species_estimate, m, num_boostrap_samples){
  n = max(m)
  for(i in 1:num_boostrap_samples){
    attribute_name = sprintf("Boostrap sample %d", i)
    n1 = length(attributes(new_species_estimate)[[attribute_name]])
    if (n1 < n){
      n = n1
    }
  }
  tmp.df = data.frame(1:n)
  for(i in 1:num_boostrap_samples){
    attribute_name = sprintf("Boostrap sample %d", i)
    tmp.df$attribute_name = attributes(new_species_estimate)[[attribute_name]][1:n]
    names(tmp.df)[names(tmp.df) == "attribute_name"] <- attribute_name
  }
  lower = apply(tmp.df, 1, function(x){quantile(x, 0.025)})
  upper = apply(tmp.df, 1, function(x){quantile(x, 0.975)})
  return(data.frame(m=1:n, lower=lower, upper=upper))
}

par(mfrow = c(2, 2))
for (n1 in names(samples_1K[c(1,2,3,4)])){
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  m =  c(freq_table %*% as.numeric(names(freq_table)))
  m = seq(1,  m * 1.3)
  
  plot(range(m), c(0, 350), type='n',
       ylab='Estimated unseen species',
       xlab='# of future samples',
       main=n1)
  lines(m, new_species_GT[[n1]], col=tab.green)
  lines(m, new_species_ET[[n1]], col=tab.blue)
  lines(m, new_species_PoiGamma[[n1]], col=tab.orange)
  
  tofill = BootstrapPredictionIntervals(new_species_GT[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.green.alpha, border = FALSE)
  tofill = BootstrapPredictionIntervals(new_species_ET[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.blue.alpha, border = FALSE)
  tofill = BootstrapPredictionIntervals(new_species_PoiGamma[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.orange.alpha, border = FALSE)
  addNewDataPoints(n1)
  abline(v=c(freq_table %*% as.numeric(names(freq_table))), lty=2, col='grey')
  legend('topleft', lty=c(1,1,1), pch=c(NA,NA,NA), c('GT','ET', 'PoiGamma'), col=c(tab.green, tab.blue, tab.orange))
}

