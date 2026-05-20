library(ggplot2)
source("FPG_estimator.R")
source("utils.R")
library(tidyverse)


host_names= c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost <- fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host


df = data.frame(matrix(data=NA, ncol=4, nrow=0))
names(df) = c("x","y", "z", "host")

for (i in 1:8){
  n1 = host_names[i]
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  par=PoissonGamma_MLE(freq_table, debug=F)
  C = freq_table[1] / dnbinom(1, par[1], par[2]/(1+par[2]))
  
  df <- rbind(df, data.frame(
    x = as.numeric(names(freq_table)),
    y = as.numeric(freq_table),
    z = 'Observed',
    host=n1))
  df <- rbind(df, data.frame(
    x = 1:max(df$x),
    y = C * dnbinom(1:max(df$x), size = par[1], prob = par[2] / (1 + par[2])),
    z = 'FPG model',
    host=n1))
}
 
ggplot() +
  geom_point(data = df[df$z=='Observed',], aes(x, y, color=z)) +
  geom_point(data = df[df$z=='FPG model',], aes(x, y), color='red', size=0.6) +
  geom_line(data = df[df$z=='FPG model',], aes(x, y, color=z)) +
  labs(x = "Phage species count", y = "Frequency") +
  theme_minimal() + scale_x_log10(limits = range(df[df$z=='Observed',]$x)) +
  scale_y_log10(limits = range(df[df$z=='Observed',]$y)) +
  facet_wrap(vars(host), ncol=2) + labs(color='')

ggsave("fit_fpg.pdf")
