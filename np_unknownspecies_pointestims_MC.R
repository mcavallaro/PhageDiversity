#' The R package  has an implementation 
#' of the Efron-Thisted estimator for 
#' the number of unseen species in the 
#' "next" sample
#' #' Install packages
#' install_pack <- list(devtools=FALSE,
#'                      variantprobs=FALSE)
#' 
#' # if (install_pack$devtools){
#' #   install.packages("remotes")}
#' if (install_pack$variantprobs){
#' remotes::install_github("c7rishi/variantprobs")}

source("NP_estimators_MC.R")
#' Load packages
# library(variantprobs) #estimators
# library(readr) #reading in csv
library(dplyr) #splitting into hosts

#' Read in species count data,
#' extract list with species count per host
#' species
fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost <- fulltable |> select(Host,
                                   `Phage Species`) |> 
                            nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host
#' count number of samples per host species
#' which ones are > 1000?
nosamples_host <- sapply(spec_byhost_l,
                         nrow)
samples_1K <- which(nosamples_host>999)

lambda <- .8# lambda=m/n

new_species_GT = list()
new_species_ET = list()


for (n1 in names(samples_1K[c(1,2,3,4)])){
  speccounts <- as.vector(unname(table(
    spec_byhost_l[[n1]])))
  freq_table<-c(table(speccounts))
  cat("\n", n1)
  cat(". n of species: ", sum(freq_table))
  cat(". n of isolates: ", c(freq_table %*% as.numeric(names(freq_table))))
  m = seq(1, sum(freq_table) * 1.3)
  new_species_GT[[n1]] = good_toulmin(freq_table, m )
  new_species_ET[[n1]] = efron_thisted(freq_table, m )
  
  # remove negatives:
  new_species_GT[[n1]]  = ifelse(new_species_GT[[n1]] < 0, 0, new_species_GT[[n1]]  )
  new_species_ET[[n1]]  = ifelse(new_species_ET[[n1]] < 0, 0, new_species_ET[[n1]]  )  
  # cat(n1,": expect ",
  #   sgt_Delta(counts = speccounts,
  #             m=sum(speccounts),
  #             t=lambda,adj = TRUE),
  #   "new species in",t*100,
  # "% more samples\n")  
}


par(mfrow = c(2, 2))    
for (n1 in names(samples_1K[c(1,2,3,4)])){
  speccounts <- as.vector(unname(table(
    spec_byhost_l[[n1]])))
  freq_table<-c(table(speccounts))
  cat("\n", n1)
  cat(" n of species: ", sum(freq_table))
  cat(". n of isolates: ", c(freq_table %*% as.numeric(names(freq_table))))
  m = seq(1, sum(freq_table) * 1.3)
  
  plot(range(m), c(0, 150), type='n',
       ylab='Estimated unseen species',
       xlab='# of future samples',
       main=n1)
  lines(m, new_species_GT[[n1]], col='green')
  lines(m, new_species_ET[[n1]], col='blue')
  abline(v=sum(freq_table), lty=2, col='grey')
  legend('topright', lty=1, c('GT','ET'), col=c('green', 'blue'))
}


# MTB has negative estimate for large t (e.g., =1),
# Have found that this may indicate saturation
# e.g. here: https://github.com/nf-core/chipseq/issues/133

