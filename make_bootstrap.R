#set.seed(14)
#par(mfrow=c(2,2))
#for (i in 1:4){
#' to get make-table-functions 
source("utils.R")
library(magrittr)
library(dplyr)

fulltable25 <- read.csv("data/3May2025_data.tsv",sep = "\t")
#h1 <- "Mycobacterium"
#h1<-"Vibrio"
h1<-"Escherichia"
specs0 <- fulltable25 |> filter(Host==h1) |> 
  select(vOTU) |> unlist()

specs <- sample(specs0,size = 250,
                replace = FALSE)

speccounts<-getSpeciesCount(specs)
freq_table<-getFrequencyTable(speccounts)
M <- extractM(speccounts)


plot(freq_table, log='x', 
     col='blue', pch=19,
     main=paste0("Species spectrum: ",
                  " host ",h1,", subsample size: ",250),
    xlab="Species abundance",ylab="# species with this abundance")

#' Make the correction sampler
e_spec_sampled <- function(sc=speccounts){
  prob_notsamp <- function(x){(1-x/sum(sc))^sum(speccounts)}
e_bt_spec <- length(sc)-sum(sapply(sc,prob_notsamp))
return(e_bt_spec)}

#plot(sapply(1:10000,f1))
specs_in_bt <- rep(-1,10)
for(i in 1:10){ # wrap FF's code to plot the frequency table from 10 bootstrap samples
#bt_naive <- sample(spec_byhost_l[[n1]]$vOTU,n)
#cat("bt # species: ",length(btn))
#' Bootstrap set
#' step 1: make augmented species set
sc_bt <- speccounts #we start with our species vector
ecount1 <- e_spec_sampled(speccounts) #no expected sampled species
run1 <- 0 #iteration counter
resamp1 <- c("weighted","uniform")[1]

while (ecount1 < length(speccounts)){
run1 <- run1 + 1
#add a new spec w. resampled abundance
if (resamp1=="uniform"){
newspec <- sample(speccounts,1)}
if (resamp1=="weighted"){
sc1 <- speccounts 
prob_notsamp <- function(x){(1-x/sum(sc1))^sum(sc1)}
prob1 <- sapply(sc1,prob_notsamp)
newspec <- sample(x = sc1,
                  1,prob = prob1)
}
names(newspec) <- paste0("votuA_",run1)
sc_bt <- c(sc_bt,newspec)
ecount1 <- e_spec_sampled(sc_bt)
}
cat("E(species lost in bt):", 
    length(speccounts)-e_spec_sampled(speccounts),"\n",
    "Added species: ",
    ecount1-e_spec_sampled(speccounts),"\n",
    "Expected sampled species BT - observed:",
    ecount1 - length(speccounts))
#' Make phage set from the augmented species set
counts2phages <- function(x){rep(names(x),x)}
phageset_augment <- NULL
for (j in seq(along=sc_bt)){
phageset_augment <- c(phageset_augment,
                      counts2phages(sc_bt[j]))}
#' Make bootstrap sample from this set
#' same size as the original data set
phageset_bt <- sample(phageset_augment,
                      size = length(specs),
                      replace = TRUE)
#' Record # species in the bt set
cat("diff species bt vs real: ",
    length(unique(phageset_bt))-length(speccounts),
    "\n")
specs_in_bt[i] <- length(unique(phageset_bt))-length(speccounts)
speccounts_bt<-getSpeciesCount(phageset_bt)
freq_table_bt<-getFrequencyTable(speccounts_bt)


#lines(freq_table_bt) # we do not plot this
}


# Let us add the freq_table from the more naive bootstrap approach:
# used to make the figure ("new_species_projection.png" in Overleaf)

for(i in 1:10){
speccounts_bt_naive <- sample(speccounts, length(speccounts), replace=T)
phages_bt <- sample(specs,length(specs),replace = TRUE)
freq_table_bt_naive<-getFrequencyTable(speccounts_bt_naive)
freq_table_bt_ind <- getFrequencyTable(getSpeciesCount(phages_bt))
lines(freq_table_bt_naive, col='red')
lines(freq_table_bt_ind)
}

points(freq_table, col='blue', cex=2)


#' add diversity in resample
#' 
#' Add 10 resamples from whole group
for(i in 1:10){
speccounts_subsamp <- getSpeciesCount(sample(specs0,
                                             size = 250,
                             replace = FALSE))
freq_table_subsamp <-getFrequencyTable(speccounts_subsamp)

lines(freq_table_subsamp, col='orange')
}
       
legend("topright",
       c("ET individual bt", "ET species bt", 
         "ET subsample",
         "ET reps"),
       pch=c(NA,NA,1,NA),
       col=c('black', 'red', 'blue',"orange"),
       lty=c(1,1,NA,1)) 
#}