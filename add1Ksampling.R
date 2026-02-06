#' this script allows to predict
#' how many new species would be in
#' an additional 1K sampled phages for
#' an "unknown" species - which we model
#' as a random of the 8 observed species


#' to get make-table-functions 
source("utils.R")
library(magrittr)
library(dplyr)
library(janitor)
#' import functions for non-parametric estimates
source("nonparam_estimators.R")
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

#' Get distribution of number of species in a "fresh"
#' 1K sample of a non-sampled host 
#' by subsampling a random present host species
#' and then a random 1K subset

get_1K <- function(size=1000){
  rhost <- sample(samples_1K,1)
  rsamp <- sample(spec_byhost_l[[rhost]]$vOTU,
                  size = size,replace = FALSE) 
  return(length(unique(rsamp)))                
    }

#' re(sub)sample nreps times
nreps <- 10000
dist_1Krand <- replicate(nreps,get_1K())  
summary(dist_1Krand)

#' Plot this, scale to percentage of
#' sims in histogram bin
library(ggplot2)
scal1 <- 100/nreps #' percentage
plotdata <- as_tibble(x=dist_1Krand)

ggplot(data=plotdata) + 
  geom_histogram(aes(x = value,
                     y = after_stat(count)*scal1),
                 binwidth = 20) +
  labs(x="# phage species in 1K samples",
       y="% of resampled sets") + 
  expand_limits(x=c(0,1000))
