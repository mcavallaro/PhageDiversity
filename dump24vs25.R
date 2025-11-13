#' Measure whether the spec difference predicted
#' vs. between 
#' 2024 and 2025 is signif.
#' 

library(magrittr)
library(dplyr)
#library(tidyverse)

#' to get make-table-functions 
source("utils.R")

# import functions for non-paramteric estimates
source("NP_estimators_MC.R")
# import functions for paramteric estimates
source("Par_estimators_MC.R")
# import functions for PYP estimates
source("pyp_EB_inference_fun.R")
#' Import validation errors
load("intval_n100_train4_allres_2025.RData")
#' import data
fulltable24 <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv",
                        check.names = FALSE)
fulltable25 <- read.csv("data/3May2025_data.tsv",sep = "\t")
fulltable25 <- fulltable25 |> mutate("in2024"=(Accession %in% fulltable24$Accession))

species_analysed <- names(valid_res) 

for (h1 in species_analysed){
data1 <- fulltable25 |> filter(Host==h1) |> select(in2024,vOTU)  
new_isos <- nrow(data1) - sum(data1$in2024) 
#' compute NAE for predicting dump May 2025 from dump Sept 2024  
speccounts<-getSpeciesCount(data1[data1$in2024,"vOTU"])
freq_table<-getFrequencyTable(speccounts)
pred_u_may25 <- efron_thisted(freq_table = freq_table,
                            m=new_isos) 
obs_u_may25 <- length(setdiff(unique(data1$vOTU),
                       unique(data1$vOTU[data1$in2024])))
obs_nae <- abs(pred_u_may25 - obs_u_may25)/obs_u_may25
#' Extract matching NAE 
valid_err <- valid_res[[h1]][rownames(valid_res[[h1]])=="NAE:",] 

#perform montecarlo test
cat(h1,": ",ecdf(c(valid_err[,"ET:0.pred"],
                   obs_nae))(obs_nae),"\n")
}
                                                       






 