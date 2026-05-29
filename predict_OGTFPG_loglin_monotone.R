library(magrittr)
library(dplyr)
library(tidyverse)
library(iNEXT) #For intra- and extrapolation
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
#nosamples_host <- sapply(spec_byhost_l, nrow)
#nospecies_host <- sapply(spec_byhost_l, function(x){nrow(unique(x))})
#samples_1K <- which(nosamples_host>999)

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

new_species_SGT <- tibble()

host_names<-c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

for (n1 in host_names ){
  # speccounts <- as.vector(unname(table(spec_byhost_l[[n1]])))
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  spec_abund <- spec_byhost_l[[n1]] %>% table() %>% sort(decreasing = T) %>% as.numeric()
  
  cat("\n", n1)
  cat(". n of species: ", sum(freq_table), " ", length(speccounts))
  m =  c(freq_table %*% as.numeric(names(freq_table)))
  cat(". n of isolates: ", m, " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
  m = seq(1, m * 1.3)

  temp1 <- SGT(freq_table, m)
  temp2 <- FisherPoissonGammaWrapper(freq_table,m)
  #' add iNEXT rarefication curve
  raref_c1 <- iNEXT(x = spec_abund,
                    q = 0,
                   endpoint = sum(speccounts) + max(m),
                   size = (m + sum(speccounts)),
                   se = 0)
  temp6 <- raref_c1$iNextEst$size_based$qD
  temp6 <- temp6[raref_c1$iNextEst$size_based$m %in% (m + sum(speccounts))]
  #' raref_c$iNextEst$size_based$qD is the
  #' allelic richness estimates (q=0
  #' Hill numbers) for the 
  #' rarefied species accumulation curve. 
  #' We compare to their extrapolation
  #' 2nd option: subsampling 100 samples for 
  #' some sizes, let's take jumps of 25
  #' which means ~ 40 to ~ 100 sizes
  #' We will focus on the first 20 classes
  #' as in 

#' 
raref_count <- function(k,specv){
  subsa <- sample(unname(unlist(specv)),
                  k,replace = FALSE)
  subsa <- getSpeciesCount(subsa)
  return(length(subsa))
}  
raref_mean <- function(k,specv,iter=100){
  mean(replicate(n = iter,
                 expr = raref_count(k,specv)))
}
raref_step <- 25
raref_c <- tibble(steps=
                   seq(raref_step,
                     sum(speccounts),
                     raref_step),
                 subsamp_specc=
                sapply(seq(raref_step,
                      sum(speccounts),
                      raref_step),
                  raref_mean,
                  specv=spec_byhost_l[[n1]])
                 )


#specacc <- specaccum(speccounts,
#                     method = "rarefaction")
#ev <- log(specacc$individuals)
#resp <- specacc$sites

lm0 <- lm(subsamp_specc ~ log(steps),data=raref_c)
lm1 <- lm(log(subsamp_specc) ~ log(steps),data=raref_c)
lm1$coefficients
sloglm_specacc <- function(x,
                           addpresent=FALSE,
                           correct1){
  t1 <- ifelse(addpresent,sum(speccounts),0)
  lm0$coefficients[2]*log(x+t1)+lm0$coefficients[1]}
loglm_specacc <- function(x,addpresent=FALSE){
  t1 <- ifelse(addpresent,sum(speccounts),0)
  exp(lm1$coefficients[2]*log(x+t1)+lm1$coefficients[1])}
#lm1$coefficients[2]*log(x+t1)+lm1$coefficients[1]}



#' assess fit of log and log-log models for
#' subsampled data
pdf(paste0("lmfits_logmodels_",n1,".pdf"))
plot(raref_c$steps,raref_c$subsamp_specc,
     main=paste(n1,": fit semi-log and log-log to SAC"),
     xlab="Number of samples",ylab="Number of species",
     pch=4,cra=c(1,1),cex=0.5)
curve(from=1,
      to=sum(speccounts),
      expr = loglm_specacc,
      add=TRUE,col="red")
curve(from=1,
      to=sum(speccounts),
      expr = sloglm_specacc,
      add=TRUE,col="blue")
dev.off()

temp3 <- loglm_specacc(m,
                       addpresent = TRUE) -
          loglm_specacc(0,
                         addpresent = TRUE) +
         sum(freq_table)
temp4 <- sloglm_specacc(m,
                        addpresent = TRUE) -
  sloglm_specacc(0,
                addpresent = TRUE) +
  sum(freq_table)
#' Add reimplementation of Ugland's approach
source("ugland_logmodel.R")
temp5 <- fit_slog_spec(freq_table = freq_table,
                       m=m,onlym = FALSE)

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
             "cm-mod. OSW" = sum(freq_table) + make_est_m(temp1),
             "loglog model" = temp3, 
             "log(sample size) model" = temp4,
             "Ugland's semilog" = temp5,
             "CJ" = temp6,
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


    
print(plt1 + facet_wrap(facet=vars(host),
                  nrow=3,ncol=3,
                  scales="free"))
ggsave(filename = "predict_new_logmodels_nobt.pdf")

