library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
# import functions for non-paramteric estimates
source("nonparam_estimators.R")

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

new_species_ET<-tibble()

host_names<-c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

for (n1 in host_names ){
  # speccounts <- as.vector(unname(table(spec_byhost_l[[n1]])))
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  cat("\n", n1)
  n_actual_species = sum(freq_table)
  cat(". n of species: ", sum(freq_table), " ", length(speccounts))
  m =  c(freq_table %*% as.numeric(names(freq_table)))
  cat(". n of isolates: ", m, " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
  m = seq(1, m * 1.3)

  temp1<-efron_thisted(freq_table, m, 
                         max_terms = 10) + n_actual_species

  make_est_m<-function(v){vm <- v
                          vm[1] <- max(0,v[1]) 
                          vm[2] <- min(max(vm[1],v[2]), 2*vm[1]) 
                          for (i in 3:length(v)){
                          vm[i] <- min(max(vm[i-1], v[i]),
                                       2*vm[i-1] - vm[i-2])  
                          }
                          return(vm)}

  t1<-tibble(m = m,
             ET = temp1,
             "mod. ET" = make_est_m(temp1),
             "host" = rep(n1, length(m)),
             n = length(spec_byhost_l[[n1]][[1]]))
  

  t1<-t1 |> pivot_longer(cols = c(-m,-host,-n),
                         names_to = "estim.")   
  t1<-t1 |> mutate(m = as.integer(m))
  t1<-t1 |> mutate(estim. = factor(estim.))

  new_species_ET<-bind_rows(new_species_ET, t1)
}

data1<-new_species_ET[new_species_ET$estim. == 'mod. ET',]

labels<-data1 %>% group_by(host) %>% slice_tail(n=1)

data1 |> ggplot() + geom_line(
    aes(x = m, y = value, colour = host)
  ) +
  geom_text_repel(
    data = labels,
    aes(x = m, y = value, label = host),
    box.padding = 0.1,
    point.padding = 0.5,
    max.overlaps = Inf,
    segment.color = 'white',
    fontface = "italic",
    size = 6
  ) +
  labs(y = "Number of phage species per host", x = "Number of future individuals observed (sampling effort)") + 
  theme_minimal() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13))


ggsave(filename = "cover.png", bg = 'white')
