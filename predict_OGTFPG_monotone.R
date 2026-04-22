library(magrittr)
library(dplyr)
library(tidyverse)
# import functions for non-parametric estimates
source("nonparam_estimators.R")
# import FPG estimator
source("FPG_estimator.R")
# import functions for bootstrap and other utils
source("utils.R")

source("import_data.R")

# import data
fulltable <- read.csv("data/3May2025_data.tsv",sep = "\t")
spec_byhost <- fulltable |> select(Host,vOTU) |> nest_by(Host)
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
# 1 vOTU_0044      
# 2 vOTU_0044
# 3 vOTU_0045
# 4 vOTU_0048
# 5 vOTU_0049
# 6 vOTU_0048
# 7 vOTU_0049
# 8 vOTU_0052
# 9 vOTU_0054
# 10 vOTU_0056
# #  1,634 more rows
# #  Use `print(n = ...)` to see more rows

# plot(sort(getSpeciesCount(spec_byhost_l[["Pseudomonas"]]), decreasing = T), log = 'xy')

# Estimate the missing species from the Sept2024 data and see how it compares with 
# the last dataset

new_species_SGT<-tibble()

host_names<-c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

for (n1 in host_names ){
  # speccounts <- as.vector(unname(table(spec_byhost_l[[n1]])))
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  cat("\n", n1)
  cat(". n of species: ", sum(freq_table), " ", length(speccounts))
  m =  c(freq_table %*% as.numeric(names(freq_table)))
  cat(". n of isolates: ", m, " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
  m = seq(1, m * 1.3)

  temp1 <- SGT(freq_table, m)
  temp2 <- FisherPoissonGammaWrapper(freq_table,m)
  
  make_est_m<-function(v){vm <- v
                          vm[1] <- max(0,v[1]) 
                          vm[2] <- min(max(vm[1],v[2]), 2*vm[1]) 
                          for (i in 3:length(v)){
                          vm[i] <- min(max(vm[i-1], v[i]),
                                       2*vm[i-1] - vm[i-2])  
                          }
                          return(vm)}

  t1<-tibble(m = m,
             FPG = sum(freq_table) + temp2,
             "cm-mod. OGT" = sum(freq_table) + make_est_m(temp1),
             "host" = rep(n1, length(m)),
             n = length(spec_byhost_l[[n1]][[1]]))
  

  t1<-t1 |> pivot_longer(cols = c(-m,-host,-n),
                         names_to = "Estimator")   
  t1<-t1 |> mutate(m = as.integer(m))
  t1<-t1 |> mutate(estim. = factor(Estimator))

  new_species_SGT <- bind_rows(new_species_SGT,t1)
}

library(ggforce)
data1 <- new_species_SGT
plt1 <- data1 |> ggplot() + geom_line(aes(x = m,y = value,
                              colour = Estimator)) +
        geom_vline(aes(xintercept = n),
                   linetype = 2
                   ) +
  labs(y="#present + #predicted species") 

pdf('predict_newspec_OGTFPG_1plotpage.pdf') 
for (i in 1:8){
  print(plt1 + facet_wrap_paginate(facet=vars(host),
                                   nrow=1,ncol=1,
                                   scales="free",
                                   page = i))
}
dev.off() 
    
print(plt1 + facet_wrap(facet=vars(host),
                  nrow=3,ncol=3,
                  scales="free"))
ggsave(filename = "predict_new_ogt_fpg_nobt.pdf")

