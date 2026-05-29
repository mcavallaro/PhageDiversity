library(magrittr)
library(dplyr)
library(tidyverse)
#library(iNEXT) #For intra- and extrapolation
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

new_species_3 <- tibble()
new_species_all <- tibble()


host_names<-c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

host_sizes <- rep(-1,8)
names(host_sizes) <- host_names

current_spec <- rep(-1,8)
names(current_spec) <- host_names

for (n1 in host_names ){
  # speccounts <- as.vector(unname(table(spec_byhost_l[[n1]])))
  speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
  freq_table<-getFrequencyTable(speccounts)
  spec_abund <- spec_byhost_l[[n1]] %>% table() %>% sort(decreasing = T) %>% as.numeric()
  
  cat("\n", n1)
  cat(". n of species: ", sum(freq_table), " ", length(speccounts))
  m =  c(freq_table %*% as.numeric(names(freq_table)))
  host_sizes[n1] <- m
  current_spec[n1] <- sum(freq_table)
  cat(". n of isolates: ", m, " ", length(spec_byhost_l[[n1]][[1]]), " ", sum(speccounts))
  m = seq(1, m * 1.3)

  temp1 <- SGT(freq_table, m)
  temp2 <- FisherPoissonGammaWrapper(freq_table,m)
  temp3 <- ChaoJost(freq_table,m)

make_est_m<-function(v){vm <- v
                          vm[1] <- max(0,v[1]) 
                          vm[2] <- min(max(vm[1],v[2]), 2*vm[1]) 
                          for (i in 3:length(v)){
                          vm[i] <- min(max(vm[i-1], v[i]),
                                       2*vm[i-1] - vm[i-2])  
                          }
                          return(vm)}

  t1<-tibble(m = m,
             FPG =  temp2,
             "OSW" =  temp1,
             "cm-mod. OSW" =  make_est_m(temp1),
             "CJ" = temp3, 
             "host" = rep(n1, length(m)),
             n = length(spec_byhost_l[[n1]][[1]]))
  #' to add Ugland's et al semi-log
  source("ugland_logmodel.R")
  temp5 <- fit_slog_spec(freq_table = freq_table,
                         m=m,onlym = TRUE)
  t2 <- bind_cols(t1,
                  "Ugland semilog"= temp5)

  t1<-t1 |> pivot_longer(cols = c(-m,-host,-n),
                         names_to = "Estimator")   
  t1<-t1 |> mutate(m = as.integer(m))
  t1<-t1 |> mutate(estim. = factor(Estimator))

  t2 <- t2 |> pivot_longer(cols = c(-m,-host,-n),
                         names_to = "Estimator")   
  t2 <- t2 |> mutate(m = as.integer(m))
  t2 <- t2 |> mutate(estim. = factor(Estimator))
  
#' Add bootstrap values for parameters
  #Make boostrap for different m
  bt_mod <- function(est=SGT){
    bt1<-sample(speccounts, length(speccounts), 
                replace=T)
    freqtab<-getFrequencyTable(bt1)
    out1<-est(freqtab, m)
    out2<-make_est_m(out1)
    return(c(out1, out2))
  }
  bt_est <- function(est){
    bt1<-sample(speccounts, length(speccounts), 
                replace=T)
    freqtab<-getFrequencyTable(bt1)
    out1<-est(freqtab, m)
    return(out1)
  }
  #' Add bt for OSW and c-m mod OSW
  bt_est_df<-replicate(100, bt_mod(), simplify = TRUE)
  bt_est_qu<-t(apply(bt_est_df, 1,
                    function(v){quantile(v,probs=c(0.025,0.975))}))
  
  bt_est_qu2<-tibble(qlow = bt_est_qu[,1],
                    qhigh = bt_est_qu[,2],
                    m2 = rep(m,2),
                    estim. = c(rep("OSW", length(m)), 
                               rep("cm-mod. OSW", length(m)))
  )
  bt_est_qu3 <- bt_est_qu2
  #' add bt for CJ
  bt_est_df<-replicate(100, bt_est(ChaoJost), 
                       simplify = TRUE)
  bt_est_qu<-t(apply(bt_est_df, 1,
                     function(v){quantile(v,probs=c(0.025,0.975))}))
  
  bt_est_qu2<-tibble(qlow = bt_est_qu[,1],
                     qhigh = bt_est_qu[,2],
                     m2 = m,
                     estim. = rep("CJ", length(m)))
  
  bt_est_qu3 <- bind_rows(bt_est_qu3,bt_est_qu2)
  
  #' add bt for FPG
  bt_est_df<-replicate(100, 
                       bt_est(FisherPoissonGammaWrapper), 
                       simplify = TRUE)
  bt_est_qu<-t(apply(bt_est_df, 1,
                     function(v){quantile(v,probs=c(0.025,0.975))}))
  
  bt_est_qu2<-tibble(qlow = bt_est_qu[,1],
                     qhigh = bt_est_qu[,2],
                     m2 = m,
                     estim. = rep("FPG", length(m)))
  
  bt_est_qu3 <- bind_rows(bt_est_qu3,bt_est_qu2)
  #' glue to non-bt parts
  bt_est_qu3 <- bt_est_qu3 |> rename(m=m2)
  t1 <- inner_join(t1, bt_est_qu3, by = c("m", "estim."))
  #' add to all host genera df
  new_species_3 <- bind_rows(new_species_3,t1)
  new_species_all <- bind_rows(new_species_all,t2) 
}

#' With bt, 3 estimators
data1 <- new_species_3
data1 |> ggplot() + geom_line(aes(x = m,y = value,
                                  colour = estim.,
                                  linetype = estim.)) +
  geom_ribbon(aes(x = m, ymin = qlow, ymax = qhigh,
                  colour = estim.,fill=estim.,
                  linetype=estim.),
              alpha = .05) +
  geom_vline(aes(xintercept = n),
             linetype = 2
  ) +
  labs(y="#predicted species") + 
  facet_wrap(facet=vars(host),
             nrow=4,ncol=2,scales="free") +
  theme_classic()

ggsave(filename = "predict_new_3estim.pdf")

#' Without bt, 3 estimators
library(ggforce)
data1 <- new_species_3
plt1 <- data1 |> ggplot() + geom_line(aes(x = m,y = value,
                              colour = Estimator)) +
        geom_vline(aes(xintercept = n),
                   linetype = 2
                   ) +
  #labs(y="#present + #predicted species") 
  labs(y="#predicted species in additional sample")

print(plt1 + facet_wrap(facet=vars(host),
                  nrow=3,ncol=3,
                  scales="free")) +
  theme(axis.text.x = element_text(size = 8))

ggsave(filename = "predict_new_3estim_nobt.pdf")

#' w/o bt, adding Ugland et al. semi-log
data1 <- new_species_all
plt1 <- data1 |> ggplot() + geom_line(aes(x = m,y = value,
                                          colour = Estimator)) +
  geom_vline(aes(xintercept = n),
             linetype = 2
  ) +
  #labs(y="#present + #predicted species") 
  labs(y="#predicted species in additional sample")

print(plt1 + facet_wrap(facet=vars(host),
                        nrow=3,ncol=3,
                        scales="free")) +
  theme(axis.text.x = element_text(size = 8))

ggsave(filename = "predict_new_3estimUgl_nobt.pdf")

#' Quantify results
for (n1 in host_names){
tempA <- data1 |> filter(host==n1,
                         m==host_sizes[n1]) |> 
  summarise(minv = min(value)) |>
  as.numeric()
cat(n1,": ",tempA
    ,", (",janitor::round_half_up(tempA/current_spec[n1],digits = 3)* 100,
    "%) \n")

}
save(new_species_SGT,file="predictnew_data.RData")
