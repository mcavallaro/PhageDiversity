library(readr) #reading in csv
library(dplyr) #splitting into hosts

#' Sept 2004 data set
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
#' remove unknown entry
samples_1K <- samples_1K[-6]
#lambda <- c(.5,.8,1,1.5) #lambda=m/n
#' Load fcts
source("pyp_EB_inference_fun.R")
source("utils.R")
#' Perform alpha,theta ML optim for PYP
mlparam_pyp <- vector("list",length(samples_1K))
names(mlparam_pyp)<- names(samples_1K)
#' double loop: compute MLE over species
#' and three optim starts

#' Take 3 starting values for a(lpha) and t(heta)
startv <- list(c("a"=0.5,"t"=1),
               c("a"=0.99,"t"=25),
               c("a"=0.7,"t"=0.01))
#' currently, issue for t negative (at least for 
#'   a=0.95,  t=-0.75) 
for (n1 in names(samples_1K)){
mlparam_pyp[[n1]] <- vector("list",3)
names(mlparam_pyp[[n1]]) <- sapply(startv,
      function(v){paste0(v,collapse=",")}) 
for (i in 1:3){
mlparam_pyp[[n1]][[i]] <- PYP_MLE(extractM(
  table(spec_byhost_l[[n1]])),
  start_alpha = startv[[i]]["a"],
  start_theta = startv[[i]]["t"])
}
#mlparam_pyp[[n1]] <- PYP_MLE(extractM(table(spec_byhost_l[[n1]])))

}
#' run expected u_nm under first params
#' 
#' PYP estimates for the exact number of additional phages 
#' sampled on each host 
#' 
source("2025data.R")
#add_samples <- 500
res1 <- sapply(1:6,function(i){
uhat_pyp(alpha = mlparam_pyp[[i]][[1]]["alpha"],
         theta = mlparam_pyp[[i]][[1]]["theta"],
         m = Summary$new_samples[Summary$host==names(samples_1K)[i]],
         M = extractM(table(spec_byhost_l[[
           names(samples_1K)[i]]])))}
)
names(res1) <- names(samples_1K)
#' GT estimates
source("NP_estimators_MC.R")
res2 <- sapply(1:6,
               function(i){good_toulmin(
                 freq_table = getFrequencyTable(getSpeciesCount(spec_byhost_l[[
                   names(samples_1K)[i]]])),
                 m = Summary$new_samples[Summary$host==names(samples_1K)[i]])})
names(res2) <- names(samples_1K)
#' Efron-Thisted
res3 <- sapply(1:6,
               function(i){efron_thisted(
                 freq_table = getFrequencyTable(getSpeciesCount(spec_byhost_l[[
                   names(samples_1K)[i]]])),
                 m = Summary$new_samples[Summary$host==names(samples_1K)[i]])})
names(res3) <- names(samples_1K)
#' lambda=m/n (fraction of new sample relative to old sample)
res_full <- rbind("pyp"=res1,"GT"=res2,"ET"=res3,
                  "obs"=Summary$new_species,
                  "lambda"=Summary |> mutate(
                    lambda=number.of.samples.2025/number.of.samples.2024-1) |> 
                    select(lambda) |> as.vector() |> unlist())
print(res_full)
