library(knitr)
library(janitor) #rounding up at 0.5
library(tidyverse)

#' Load data, bring into 
#' tibble form, only second half 
load("intval_n500_train5_rawdist_2025.RData")
res25 <- bind_cols(valid_res[[1]][501:1000,],
            host=rep(names(valid_res)[1],
            nrow(valid_res[[1]])/2))
for (i in 2:length(valid_res)){
  res25 <- bind_rows(res25,
              bind_cols(valid_res[[i]][501:1000,],
                        host=rep(names(valid_res)[i],
                        nrow(valid_res[[i]])/2)))
}
#' Only retain the NAEs
load("intval_n500_train5_rawdist_2024.RData")
res24 <- bind_cols(valid_res[[1]][501:1000,],
          host=rep(names(valid_res)[1],
                   nrow(valid_res[[1]])/2))
for (i in 2:length(valid_res)){
  res24 <- bind_rows(res24,
                     bind_cols(valid_res[[i]][501:1000,],
                               host=rep(names(valid_res)[i],
                                      nrow(valid_res[[i]])/2)))
}


#' Prepare for plots
plotout <- FALSE#set whether to plot outliers
bp_suffix <- ifelse(plotout,"","_noout")
#' 24 
plotd <- res24 |> pivot_longer(cols = !host,names_to = "stat",
                               values_to = "value")
plotd <- plotd |> separate_wider_delim(cols = stat,delim = ":",
                                       names = c("stat","training"))
plotd <- plotd |>  mutate(stat = fct_relevel(stat, 
                                             c("GT","ET","PYP","FPG")))
#' Check: which values are NA
plotd |> filter(is.na(value)) |> 
  group_by(host,stat,training) |> count()
#' issue: Strep 0.25 is nearly all missing
#' Check: which values are >2
plotd |> filter(value>2) |> 
  group_by(host,stat,training) |> count()
#' result: only GT for 0.25 and 3.5 and one 
#' value for FPG (which we remove) 
#' plot: facets=hosts, no GT for 0.25 and 3.5
#' remove the single > 2 value
ggplot(data = plotd |> 
         filter(!(training %in% c("0.25","0.35") &
                  stat=="GT")) |> 
         filter(value<=2),
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot(outlier.size = .4,outliers=plotout) + 
  labs(title = 2024,
                        x="training set size",
                        y="NAE") + 
  facet_wrap(vars(host),
             nrow = 4,
             scales="free_y")
ggsave(paste0("intval_2024_byhost",bp_suffix,".pdf"),
       height = 6,width = 5)

#' same plot, by stats, across training sizes
ggplot(data = plotd |> 
         filter(!(training %in% c("0.25","0.35") &
                    stat=="GT")) |> 
         filter(value<=2),
       aes(x=host,y=value,fill=training)) + 
  geom_boxplot(outlier.size = .4,outliers=plotout) + 
  labs(title = 2024,
       x="Host",
       y="NAE") + 
  theme(axis.text.x = element_text(angle = 90,
        vjust = 0.5, hjust=1)) +
  facet_wrap(vars(stat),
             nrow = 2,
             scales="free_y")
ggsave(paste0("intval_2024_bystat",bp_suffix,".pdf"),
       height = 4,width = 8)

plotd_24 <- plotd
#2025

#' Prepare for plots
plotd <- res25 |> pivot_longer(cols = !host,names_to = "stat",
                               values_to = "value")
plotd <- plotd |> separate_wider_delim(cols = stat,delim = ":",
                                       names = c("stat","training"))
plotd <- plotd |>  mutate(stat = fct_relevel(stat, 
                                             c("GT","ET","PYP","FPG")))
plotd <- plotd |>  filter(training!="0.pred")

#' Check: which values are NA
plotd |> filter(is.na(value)) |> 
  group_by(host,stat,training) |> count()
#' issue: Strep 0.25 is all missing
#' Check: which values are >2
plotd |> filter(value>2) |> 
  group_by(host,stat,training) |> count()
#' result: only GT for 0.25 and 3.5 
#' plot: facets=hosts, no GT for 0.25 and 3.5
ggplot(data = plotd |> 
         filter(!(training %in% c("0.25","0.35") &
                    stat=="GT")),
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot(outlier.size = .4,outliers = plotout) + 
  labs(title = 2025,
       x="training set size",
       y="NAE") + 
  facet_wrap(vars(host),
             nrow = 4,
             scales="free_y")
ggsave(paste0("intval_2025_byhost",bp_suffix,".pdf"),
       height = 6,width = 5)

#' same plot, by stats, across training sizes
ggplot(data = plotd |> 
         filter(!(training %in% c("0.25","0.35") &
                    stat=="GT")),
       aes(x=host,y=value,fill=training)) + 
  geom_boxplot(outlier.size = .4,outliers=plotout) + 
  labs(title = 2025,
       x="Host",
       y="NAE") + 
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1)) +
  facet_wrap(vars(stat),
             nrow = 2,
             scales="free_y")
ggsave(paste0("intval_2025_bystat",bp_suffix,".pdf"),
       height = 4,width = 8)


plotd_25 <- plotd
save(plotd_25,plotd_24,file="plot_data.RData")
