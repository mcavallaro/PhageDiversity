library(knitr)
library(janitor) #rounding up at 0.5
library(tidyverse)


#' Load data, bring into 
#' tibble form, only second half 
load("intval_n500_fixedm_2025.RData")
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

#' Prepare for plots
plotd <- res25 |> pivot_longer(cols = !host,names_to = "stat",
                               values_to = "value")
plotd <- plotd |> separate_wider_delim(cols = stat,delim = ":",
                                       names = c("stat","training"))

plotd <- plotd |> separate_wider_delim(cols = training,delim = "_",
                                       names = c("# training","# predicted"))

ps1 <- vector("list",length=length(valid_res))
names(ps1) <- names(valid_res)
for (n1 in names(ps1)){
ps1[[n1]] <- ggplot(data = plotd |> filter(host==n1),
       aes(x=`# predicted`,y=value,fill=stat)) + 
  geom_boxplot(outlier.size = .4) + 
  labs(title = paste("Predict new species in sample of size m=... , host",n1),
       y="NAE") + 
  facet_wrap(vars(`# training`),
             scales="free_y")}
pdf("plots_fixedm.pdf")
for (i in seq(along=ps1)){
  print(ps1[[i]])
}
dev.off()

#' w/o outliers
ps1 <- vector("list",length=length(valid_res))
names(ps1) <- names(valid_res)
for (n1 in names(ps1)){
  ps1[[n1]] <- ggplot(data = plotd |> filter(host==n1),
                      aes(x=`# predicted`,y=value,fill=stat)) + 
    geom_boxplot(outlier.size = .4,outliers = FALSE) + 
    labs(title = paste("Predict new species in sample of size m=... , host",n1),
         y="NAE") + 
    facet_wrap(vars(`# training`),
               scales="free_y")}
pdf("plots_fixedm_noout.pdf")
for (i in seq(along=ps1)){
  print(ps1[[i]])
}
dev.off()


