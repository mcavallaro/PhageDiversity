#' Massimo's implementation
source("NP_estimators_MC.R")
#' variantprobs implementation
library(variantprobs)
#' Fabian's quick & dirty implementation
#' inp can be either 
#' a frequency table (named vec, names are freq classes and
#' entries are how often that freq class appears, so number of obs. species 
#' with that count) or
#' the counts per species
#' type is string: "counts" or "freqT" 
lazy_gte <- function(inp,t,type="counts"){
  if (type=="counts"){
  out1 <- -sum((-t)^inp)}
  if (type=="freqT"){
  out1 <- -sum((-t)^as.integer(names(inp))*inp)}   
  return(out1)}

#
library(dplyr) #for splitting into hosts

#' Read in species count data,
#' extract list with species count per host
#' species
fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost <- fulltable |> select(Host,
                                   `Phage Species`) |> 
  nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host

#' Use Mycobacterium as example
spec_table <- as.vector(table(spec_byhost_l$Mycobacterium))
freq_table <- c(table(spec_table))
#' As we have simple upper and lower bounds, we assume to sample
#' as many new phages as we already have
samp_already <- sum(spec_table)
samp_new <- samp_already
#head(c(freq_table))
#sgt_Delta(counts = spec_table,t = 1/sum(spec_table),m = sum(spec_table),adj = FALSE)
test1 <- c(sgt=sgt_Delta(counts = spec_table,t = samp_new/samp_already,
                         m = samp_already,
                         adj = FALSE),
           MCimp=good_toulmin(freq_table,m = samp_already),
           FFimp_counts=lazy_gte(spec_table,
                                 t = samp_new/samp_already,
                                 type = "counts"),
           FFimp_freqT=lazy_gte(freq_table,
                                t = samp_new/samp_already,
                                type = "freqT")
           )
#' Quick check: G-T for t=1 (lambda) is in [-# sampled ind,#sampled individuals] 
#' (very broad bounds, actually #species suffices)
bounds_gte_meqn <- c(-sum(freq_table),sum(freq_table)) 
test1 >= bounds_gte_meqn[1]
test1 <= bounds_gte_meqn[2]
test1
#' Use MC's test set
#Example usage
# Suppose we have observed frequencies: 118 species seen 1 time, 74 species seen 2 times, etc.
freq_table <- c(118,	74,	44,	24,	29,	22,	20,	19,	20,	15,	12,	14,	6, 1, 6)
# Corbet -FIsher butterfly data

#' Get # species and sample size
names(freq_table) = as.character(1:length(freq_table))
cat("number of species already observed:",sum(freq_table),"\n")
samp_already <- sum(as.integer(names(freq_table))*freq_table)
cat("Number of individuals sampled already: ",samp_already,"\n")

#' compute all estimators, again for t=1
test2 <- c(sgt=sgt_Delta(r = as.integer(names(freq_table)),
                         N_r = freq_table,
                         t = 1,
                         m = sum(freq),
                         adj = FALSE),
           MCimp=good_toulmin(freq_table,m = samp_already),
           FFimp_freqT=lazy_gte(freq_table,t = 1,type = "freqT")
)

bounds_gte_meqn <- c(-sum(freq_table),sum(freq_table)) 
test2 >= bounds_gte_meqn[1]
test2 <= bounds_gte_meqn[2]
test2

