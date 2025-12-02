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
plotd <- res24 |> pivot_longer(cols = !host,names_to = "stat",
                               values_to = "value")
plotd <- plotd |> separate_wider_delim(cols = stat,delim = ":",
                                       names = c("stat","training"))
plotd <- plotd |>  mutate(stat = fct_relevel(stat, 
                                             c("GT","ET","PYP","FPG")))

#' Doublecheck: only >2 for GT
#plotd |> filter(value>2) |> select(stat) |> table()
#' Set values to at most 2
#plotd$value <- sapply(plotd$value,function(x){min(2,x)})

#' As GT is bad for lambda>1 we make 2 plots - 
#' one for lamba <=1 with all statistics and one for all
#' lambda for all but GT

#' plot for small lambda
ggplot(data = plotd |> filter(!training %in% c("0.25","0.35")),
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2024,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2024_gt.pdf")

#' plot for lambda > 1, no GT
ggplot(data = plotd |> 
         filter(!stat == "GT",
                training %in% c("0.25","0.35")),
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2024,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2024_nogt_smalltrain.pdf")


#2025

#' Prepare for plots
plotd <- res25 |> pivot_longer(cols = !host,names_to = "stat",
                               values_to = "value")
plotd <- plotd |> separate_wider_delim(cols = stat,delim = ":",
                                       names = c("stat","training"))
plotd <- plotd |>  mutate(stat = fct_relevel(stat, 
                                             c("GT","ET","PYP","FPG")))
plotd <- plotd |>  filter(training!="0.pred")

#' Doublecheck: only >2 for GT
#plotd |> filter(value>2) |> select(stat) |> table()

#' Set values to at most 2
#plotd$value <- sapply(plotd$value,function(x){min(2,x)})

ggplot(data = plotd |>
         filter(!training %in% c("0.25","0.35")),
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2025_gt.pdf")

#' plot for lambda > 1, no GT
ggplot(data = plotd |> 
         filter(!stat == "GT",
                training %in% c("0.25","0.35")),
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2025_nogt_smalltrain.pdf")


ggplot(data = plotd |> filter(stat!="GT"),
       aes(x=host,y=value,fill=training)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + theme() + 
  facet_wrap(vars(stat),
             scales="free_y") +
theme(axis.text.x = element_text(angle = 90))
ggsave("intval_2025_bytraining.pdf")



ggplot(data = plotd,
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + 
  lims(y=c(0,0.2)) +
  facet_wrap(vars(host),scales="free_y")
#ggsave("intval_2025_zoom.pdf")


