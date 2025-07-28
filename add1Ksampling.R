#' to get make-table-functions 
source("utils.R")
library(magrittr)
library(dplyr)
library(janitor)
# import functions for non-paramteric estimates
source("NP_estimators_MC.R")
#' Get data
fulltable <- read.csv("data/3May2025_data.tsv",sep = "\t")
#spec_byhost <- fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
spec_byhost <- fulltable |> select(Host,vOTU) |> nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host
#' count number of samples per host species
#' which ones are > 1000?
nosamples_host <- sapply(spec_byhost_l, nrow)
nospecies_host <- sapply(spec_byhost_l, function(x){nrow(unique(x))})
samples_1K <- which(nosamples_host>999)
samples_1K <- samples_1K[-which(names(samples_1K)=="Unspecified")]
#' Get sample sizes and number of present species
#' We can always run ET/GT, we have at least 1K 
library(tibble)
obs_plus1K_may25 <- data.frame(cbind(phage_sampsize=nosamples_host[samples_1K],
                   species_sampsize=nospecies_host[samples_1K]))
obs_plus1K_may25 <- rownames_to_column(obs_plus1K_may25,var = "host")

pred1K <- function(host,m=1000){
  n1 <- samples_1K[host]
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  unname(efron_thisted(freq_table = freq_table,
                m = m))
}

obs_plus1K_may25 <- tibble(obs_plus1K_may25,
                           pred1K=sapply(obs_plus1K_may25$host,pred1K)) 
obs_plus1K_may25 <- obs_plus1K_may25 |> mutate(perc_new = 100*pred1K/species_sampsize)
knitr::kable(obs_plus1K_may25,format = "rst",digits = 2,align = "l")