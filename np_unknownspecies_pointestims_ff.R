#' The R package  has an implementation 
#' of the Efron-Thisted estimator for 
#' the number of unseen species in the 
#' "next" sample

#' Install packages
install_pack <- list(devtools=FALSE,
                     variantprobs=FALSE)

if (install_pack$devtools){
  install.packages("remotes")}
if (install_pack$variantprobs){
remotes::install_github("c7rishi/variantprobs")}
#' Load packages

library(variantprobs) #estimators
library(readr) #reading in csv
library(dplyr) #splitting into hosts

#' Read in species count data,
#' extract list with species count per host
#' species
fulltable <- read_csv("data/phagesspeciescounts_perhostspec_Sept2024.csv")
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
for (n1 in names(samples_1K)){
speccounts <-   as.vector(unname(table(
    spec_byhost_l[[n1]])))
cat(n1,": expect ",
    sgt_Delta(counts = speccounts,
              m=sum(speccounts),
              t=lambda,adj = TRUE),
    "new species in",lambda*100,
  "% more samples\n")  
}
# MTB has negative estimate for large t (e.g., =1),
# Have found that this may indicate saturation
# e.g. here: https://github.com/nf-core/chipseq/issues/133

