# Estimate the missing species from DB24 and see how it compares with DB25

library(magrittr)
library(dplyr)

# import functions for non-paramteric estimates
source("nonparam_estimators.R")

# import functions for paramteric estimates
source("FPG_estimator.R")

# import functions for PYP estimates
source("PYP_estimator.R")

# import functions for bootstrap and other utils
source("utils.R")

source("import_data.R")

fulltable<-read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost<-fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
spec_byhost_l<-as.list(spec_byhost$data)
names(spec_byhost_l)<-spec_byhost$Host
nosamples_host<-sapply(spec_byhost_l, nrow)
nospecies_host<-sapply(spec_byhost_l, function(x){nrow(unique(x))})
samples_1K<-which(nosamples_host > 999)


new_species_GT<-list()
new_species_ET<-list()
new_species_PoiGamma<-list()
new_species_PYP<-list()

num_boostrap_samples<-100

host_names= c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

for (n1 in host_names ){
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  cat("\n", n1)
  cat(". n of species: ", sum(freq_table), " ", length(speccounts))
  m<-c(freq_table %*% as.numeric(names(freq_table)))
  cat(". n of isolates: ", m, " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
  m<-seq(1, m * 1.3)
  new_species_GT[[n1]]<-good_toulmin(freq_table, m)
  new_species_ET[[n1]]<-efron_thisted(freq_table, m, 20)
  new_species_PoiGamma[[n1]]<-FisherPoissonGammaWrapper(freq_table, m)
  M<-extractM(speccounts)
  new_species_PYP[[n1]]<-BalocchiPYPWrapper(M, m)
  
  tmp=bootstrapSpeciesCount(speccounts, num_boostrap_samples)
  #tmp=bootstrapSpeciesCountFF(speccounts, num_boostrap_samples)
  for(i in 1:num_boostrap_samples){
    sample<-tmp[[i]]
    freq_table<-getFrequencyTable(sample)
    
    attribute_name<-sprintf("Boostrap sample %d", i)
    attr(new_species_GT[[n1]], attribute_name)<-good_toulmin(freq_table, m)
    attr(new_species_ET[[n1]], attribute_name)<-efron_thisted(freq_table, m, 20)
    attr(new_species_PoiGamma[[n1]], attribute_name)<-FisherPoissonGammaWrapper(freq_table, m)
    attr(new_species_PYP[[n1]], attribute_name)<-BalocchiPYPWrapper(M, m)
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

png(filename = "new_species_projection_DB24.png", width = 1400 * 1.4, height = 1600 * 1.4, res = 300)
par(mfrow = c(4, 2),
    mar = c(5, 5, 4, 2) - 1,
    oma = c(1, 2, 0, 0))
for (n1 in host_names){
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  m<-c(freq_table %*% as.numeric(names(freq_table)))
  m<-seq(1,  m * 1.3)
  
  plot(c(0, max(m)-100), c(0, 280), type='n',
       ylab='Est. unseen species',
       xlab='Num. of future samples',
       main=n1)
  lines(m, new_species_GT[[n1]], col=tab.green)
  lines(m, new_species_ET[[n1]], col=tab.blue)
  lines(m, new_species_PoiGamma[[n1]], col=tab.orange)
  lines(m, new_species_PYP[[n1]], col=tab.purple)
  
  tofill<-BootstrapPredictionIntervals(new_species_GT[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.green.alpha, border = FALSE)
  tofill<-BootstrapPredictionIntervals(new_species_ET[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.blue.alpha, border = FALSE)
  tofill<-BootstrapPredictionIntervals(new_species_PoiGamma[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.orange.alpha, border = FALSE)
  tofill<-BootstrapPredictionIntervals(new_species_PYP[[n1]], m, num_boostrap_samples)
  polygon(c(tofill$m, rev(tofill$m)), c(tofill$lower, rev(tofill$upper)), col=tab.purple.alpha, border = FALSE)
  addNewDataPoints(n1)
  abline(v=c(freq_table %*% as.numeric(names(freq_table))), lty=2, col='grey')
  if (n1=="Escherichia"){
    legend('topright', lty=c(1,1,1,1), pch=c(NA,NA,NA,NA), c('GT','ET', 'FPG', 'PY'),
           col=c(tab.green, tab.blue, tab.orange, tab.purple), cex=0.8)
  }
}
dev.off()

