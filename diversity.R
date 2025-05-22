
library(iNEXT)
library(magrittr)
library(dplyr) #splitting into hosts

#' Read in species count data,
#' extract list with species count per host
#' species
fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
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



results = list()
for (n1 in names(samples_1K[c(1,2,3,4)])){
  species_aboundance = spec_byhost_l[[n1]] %>% table() %>% sort(decreasing = T) %>% as.numeric()
  # # n. of species  
  # species_aboundance %>% length()
  # #  of individuals 
  # species_aboundance %>% sum
  results[[n1]] = species_aboundance
}


#  this iNEXT function calculated the so-called Hill numbers (which are generalised entropies)
#  with order parameters q (q=1 corresponds to the Shannon entropy)
A = iNEXT(results, q=c(0,1,2))


# iNEXT package contains nice custom plotting functions
plot(A)

# we can tell that hosts Escherichia and Klesbiella have higher phage diversity than Mycobacterium and Pseudomonas
