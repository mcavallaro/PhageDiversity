library(knitr)
library(janitor) #rounding up at 0.5
library(tidyverse)

#' Load data, bring into 
#' tibble form, only second 100 
load("intval_n100_train4_allres_2025.RData")
res25 <- bind_cols(valid_res[[1]][101:200,],
            host=rep(names(valid_res)[1],
            nrow(valid_res[[1]])/2))
for (i in 2:length(valid_res)){
  res25 <- bind_rows(res25,
              bind_cols(valid_res[[i]][101:200,],
                        host=rep(names(valid_res)[i],
                        nrow(valid_res[[i]])/2)))
}
#' Only retain 100 NAEs
load("intval_n100_train4_allres_2024.RData")
res24 <- bind_cols(valid_res[[1]][101:200,],
          host=rep(names(valid_res)[1],
                   nrow(valid_res[[1]])/2))
for (i in 2:length(valid_res)){
  res24 <- bind_rows(res24,
                     bind_cols(valid_res[[i]][101:200,],
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
plotd |> filter(value>2) |> select(stat) |> table()
#' Set values to at most 2
plotd$value <- sapply(plotd$value,function(x){min(2,x)})

ggplot(data = plotd,
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2024,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2024.pdf")

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
plotd |> filter(value>2) |> select(stat) |> table()

#' Set values to at most 2
plotd$value <- sapply(plotd$value,function(x){min(2,x)})

ggplot(data = plotd,
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2025.pdf")

ggplot(data = plotd,
       aes(x=host,y=value,fill=training)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + 
  facet_wrap(vars(stat),scales="free_y")
ggsave("intval_2025_reord.pdf")

ggplot(data = plotd,
       aes(x=training,y=value,fill=stat)) + 
  geom_boxplot() + labs(title = 2025,
                        x="% training",
                        y="NAE") + 
  lims(y=c(0,0.2)) +
  facet_wrap(vars(host),scales="free_y")
ggsave("intval_2025_zoom.pdf")


