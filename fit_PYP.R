#' Script for assessing goodness-of-fit
#' of PYP estimator
#' We ML estimate the PYP parameters
#' and then simulate the PYP with these
library(ggplot2)
source("PYP_estimator.R")
source("utils.R")
library(tidyverse)

#' We need a PYP simulator, 
#' which we get via sequential
#' generalized Chinese Restaurant process
#' construction
simulate_2pCRP <- function(ncust,alpha,
                           theta,frqt=TRUE){
  tables <- 1 #first table has first customer
  for (i in 2:ncust){
    #sequential CRP table occupation probas
    probs <- sapply(tables,function(x){(x-alpha)/(i-1+theta)})
    probs <- c(probs,1-sum(probs))
    cust_table <- sample.int(length(tables)+1,
                             size = 1,prob = probs)
    #following either opens new table or adds customer
    if (cust_table>length(tables)){
      tables <- c(tables,1)} else {
        tables[cust_table] <- tables[cust_table] +1}
  }  
  out1 <- tables
  if (frqt){out1 <- getFrequencyTable(out1)}
  return(out1)
}  

#'
host_names= c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost <- fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host

capped1 <- 24 
df <- tibble()



library(future.apply)
plan(multisession,workers=3)
for (i in 1:8){
  n1 <- host_names[i]
  speccounts <- getSpeciesCount(spec_byhost_l[[n1]])
  freq_table <- getFrequencyTable(speccounts)
  par <- PYP_MLE(M = extractM(speccounts),
                 debug = FALSE)

cappf <- function(tab1,capped=24){
  all_counts <- sum(tab1)
  out1 <- unname(tab1[as.character(1:capped)])
  out1[is.na(out1)] <- 0
  names(out1) <- paste0("f",1:capped)
  out1 <- c(out1,"proxy"=all_counts-sum(out1))
  names(out1)[capped+1] <- paste0("f",capped+1,"+")  
  return(out1)
}
#######
  f_run <- function(freq=FALSE,capped=24){ 
    out1 <- simulate_2pCRP(sum(speccounts),
                           theta = par["theta"], 
                           alpha = par["alpha"],
                           )
    out1 <- cappf(out1,capped = capped) 
    if (freq){out1 <- out1/sum(out1)}
    return(out1)
  }

  C1 <- future_replicate(expr = f_run(),
                         n = 500)   
  C <- t(C1) #add. C1 is for code development stage
  C <- bind_cols(C,
              host=rep(n1,nrow(C)),
              type="from fitted")
  obs1 <- as_tibble(t(cappf(freq_table)))
  obs1 <-  tibble(obs1,
                 host=n1,
                 type="observed")
  C <- bind_rows(as_tibble(C),obs1)
  df <- bind_rows(df,C)
}
df1 <- df #for debugging
save(df,file="fit_PYP.RData")
#df <- 
df |> pivot_longer(f1:`f25+`,
                   names_to = "count class",
                   values_to = "count") 
df$`count class` <- factor(df$`count class`,
                           levels=c(paste0("f",1:24),"f25+"))

df2 <- 
  df |> group_by(host,type,`count class`) |> 
  summarise(med=median(count))

ggplot(data=df2) +
  geom_point(aes(x=`count class`,y=med,
                 colour=type,shape=type)) + 
  facet_wrap(vars(host),nrow=4,ncol=2) + 
  scale_y_sqrt() +
  labs(x = "Phage species count", y = "Frequency") +
  theme_minimal() 

ggsave(filename = "fit_pyp.pdf")




