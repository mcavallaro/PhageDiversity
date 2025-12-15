#' We internally validate all estimators, by taking 
#' a random subset of the database 
#' of some size (50%,65%,80%) and predicting
#' the rest of the database  

#' to get make-table-functions 
source("utils.R")

#' From a species set, extract two subsets 
#' (training and validation) of sizes size and m
#' directly from species freq table or 
#' from the Mr,n (see utils.R for details)
#' and return the training set and the # species present 
#' in the validation data set but not in the training data set
#' M= named freq table or Mr,n vector, so a table
#' showing how many species appear 1,2,...,n times in 
#' n sampled phages 
#' size: size of training subset 
#' m: size of validation set
#' type=freq_table or Mrn
#' We do this by extending the table to a vector representing
#' the species, split this randomly and then   
#' re-table both subsets
make_train_val_sets2 <- function(M,size,m=50,type="freq_table"){
if (type=="freq_table"){
  #' make a vector of characters for each species occured 
  #' that denote the number of phages found from this species 
  sampfrom <- unlist(sapply(seq(along=M),
                     function(i){rep(names(M)[i],
                                     M[i])})) 
  #' from this, make a vector of species s1,s2,... 
  #' that fit to these counts
  spec_names <- paste0("s",seq(along=sampfrom))
  sampfrom <- unlist(sapply(seq(along=sampfrom),
                     function(i){rep(spec_names[i],
                                      as.integer(sampfrom[i]))})) 
  #' reorder this dataset
  subsample <- sample(sampfrom) 
  #' split, record species in training and validation set
  train1 <-  subsample[1:size]
  specs_seen <- list("train"=unique(train1),
                     "val"=unique(subsample[-(1:size)][1:m]))
  #' split into training and validation subsets
  return(list("trainset"=
         getFrequencyTable(getSpeciesCount(train1)),
         "val_n_speccount"=length(setdiff(specs_seen$val,
                                       specs_seen$train))))
}
if (type=="Mrn"){
  #' make a vector of characters for each species occured 
  #' that denote the number of phages found from this species 
  sampfrom <- unlist(sapply(seq(along=M),function(i){rep(i,M[i])}))
  #' from this, make a vector of species s1,s2,... 
  #' that fit to these counts
  spec_names <- paste0("s",seq(along=sampfrom))
  sampfrom <- unlist(sapply(seq(along=sampfrom),
                     function(i){rep(spec_names[i],
                                      sampfrom[i])})) 
  #' reorder the dataset
  subsample <- sample(sampfrom) 
  #' split, record species in training and validation set
  train1 <-  subsample[1:size]
  specs_seen <- list("train"=unique(train1),
                     "val"=unique(subsample[-(1:size)][1:m]))
  #' record new species in val set, training set
  return(list("trainset"=extractM(getSpeciesCount(train1)),
              "val_n_speccount"=length(setdiff(specs_seen$val,
                                            specs_seen$train))))
}
}  

#use summary with NA's always recorded
summaryna <- function(v){
  out1 <- summary(v)
  if(!any(is.na(v))){
    out1 <- c(out1,"NA's"=0)
  }
  return(out1)
}

#' Use internal validation: split a dataset 
#' n times into a test and validation subset,
#' use a estimation function to predict the number of species 
#' in the validation set
#' estimation function must have argument structure and order estim(M,m) 
#' [all our wrappers have this]
#' argument estim: function for estimation
#' M= named freq table or Mr,n vector (the data) 
#' trainsize: size of training subset
#' m: size of validation set
#' type=freq_table or Mrn 
#' out: 3 options 
#' "rawdist" outputs 2-col matrix w. 
#' absolute errors (col 1) and NMAEs (col2)
#' for the n assessments, 
#' "summdist" only a summary
#' using summaryna above for each column
#' "raw": true values and estimated values   
intval2 <- function(estim1=BalocchiPYPWrapper,
                   trainsize,M,m=50,
                   type="Mrn",n=10,
                   out="raw"){
res <- cbind("abs_err"=rep(-1,n),
             "NAE"=rep(-1,n))
raw1 <- cbind("true"=rep(-1,n),
              "estim"=rep(-1,n))
#' counts number of phages  
totalsize <- switch(type,
                    "freq_table"=sum(as.integer(names(M))*M),
                    "Mrn"=sum(seq(along=M)*M))
for (i in 1:n){
split1 <- make_train_val_sets2(M,size = trainsize,m = m,
                              type=type)

estim_val <- estim1(split1$trainset,m)
#' record # new species estimate with real value
if (out=="raw"){
  raw1[i,2] <- estim_val
  raw1[i,1] <- split1$val_n_speccount  
} else {
res[i,1] <- abs(estim_val-split1$val_n_speccount)
res[i,2] <- abs(estim_val-split1$val_n_speccount)/split1$val_n_speccount}
}
out1 <- switch(out,
               "rawdist"=res,
               "summdist"=apply(res,2,summaryna),
               "raw"=raw1)
return(out1)
}
####################
####################
#' Example analysis
library(magrittr)
library(dplyr)
# import functions for non-parametric estimates
source("NP_estimators_MC.R")
# import functions for parametric estimates
source("Par_estimators_MC.R")
# import functions for PYP estimates
source("pyp_EB_inference_fun.R")
# import functions for bootstrap and other utils
#source("utils.R")
#source("2025data.R")

#' import data
fulltable24 <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv",
                      check.names = FALSE)
fulltable25 <- read.csv("data/3May2025_data.tsv",sep = "\t")
fulltable25 <- fulltable25 |> mutate("in2024"=(Accession %in% fulltable24$Accession))

for (h1 in c("Escherichia","Klebsiella","Mycobacterium","Pseudomonas",
             "Salmonella","Staphylococcus","Streptococcus","Vibrio")){
cat(h1,"--","2025 filtered: ",  
sum(fulltable25 |> filter(Host==h1)|> select(in2024)),
"vs 2024: ",nrow(fulltable24 |> filter(Host==h1)),"\n")
}
#Perform internal validations
yr <- 2025
#' we are running the analysis for the 2025 table
#spec_byhost <- fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
  if (yr==2024){
    spec_byhost <- fulltable25 |> filter(in2024==TRUE) |> select(Host,vOTU) |> nest_by(Host)}
  if (yr==2025){
    spec_byhost <- fulltable25 |> select(Host,vOTU) |> nest_by(Host)}
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host
#' count number of samples per host species
#' which ones are > 1000?
nosamples_host <- sapply(spec_byhost_l, nrow)
nospecies_host <- sapply(spec_byhost_l, function(x){nrow(unique(x))})
samples_1K <- c("Escherichia","Klebsiella",
                "Mycobacterium","Pseudomonas",
                "Salmonella","Staphylococcus",
                "Streptococcus","Vibrio")

#' run validation
#' 
valid_res <- vector("list",length = length(samples_1K))
names(valid_res) <- samples_1K

valreps <- 500#50

for (i in seq(along=samples_1K)){
n1 <- samples_1K[i]
 for (m in seq(100,500,100)){
speccounts <- getSpeciesCount(spec_byhost_l[[n1]])
freq_table <- getFrequencyTable(speccounts)
M <- extractM(speccounts)
#' approx. 80%,... for training
trainsize <- c(200,300,400,500)
#' add actual size from 2024 as training size


#' ET (safety and sanity)
temp1 <- sapply(trainsize,function(s1){intval2(estim1=efron_thisted,m=m,
         trainsize = s1,
         M = freq_table,
         type = "freq_table",n = valreps,out = "rawdist")},
  simplify = "matrix")
colnames(temp1) <- paste0("ET:",trainsize,"_",m)
valid_res[[i]] <- cbind(valid_res[[i]],temp1)
#' FisherPoissonGamma
temp1 <- sapply(trainsize,function(s1){
  intval2(estim1=FisherPoissonGammaWrapper,
         trainsize = s1,
         M = freq_table,type = "freq_table",
         n = valreps,out = "rawdist")})
colnames(temp1) <- paste0("FPG:",trainsize,"_",m)
valid_res[[i]] <- cbind(valid_res[[i]],temp1)
}}
save(valid_res,file = paste0("intval_n",valreps,"_fixedm_",yr,".RData"))
