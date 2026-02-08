#' This script computes various
#' summary tables from the output data
#' from script plots_intval_fixedm.R
#' and perform Wilcoxon tests
#' between the internal validation errors 
#' of estimators ET and FPG

library(tidyverse)
load("plot_data.RData")

#' extract estimator with lowest mean errors 
testf <- function(v,w){v[which(w==min(w))]}
meanerr <- plotd_24 |> 
  group_by(host,training,stat) |> 
  summarise(sum1=mean(value)) 
t1 <- meanerr |> reframe(bestmod=testf(stat,sum1))
cat("2024")
#tables show which estimator had lowest min error
table(t1$bestmod) 
tapply(t1$bestmod,t1$training,table)

meanerr <- plotd_25 |> filter() |> 
  group_by(host,training,stat) |> 
  summarise(sum1=mean(value)) 
t1 <- meanerr |> reframe(bestmod=testf(stat,sum1))
cat("2025")
#tables show which estimator had lowest min error
table(t1$bestmod)
tapply(t1$bestmod,t1$training,table)

#' Mean errors via training/validation ratio and estimator
plotd_24 |> group_by(training,stat) |> summarise(mean(value)) |>
  View()
plotd_25 |> group_by(training,stat) |> summarise(mean(value)) |>
  View()
#' only for ET
plotd_25 |> group_by(training,stat) |> 
  summarise(m= mean(value)) |> filter(stat=="ET") |> select(training,m)

#' Compute Wilcoxon tests
wilcx_p <- matrix(-1,nrow=length(unique(plotd_25$host)),
                     ncol=length(unique(plotd_25$training)),
                     dimnames=list(unique(plotd_25$host),
                                   unique(plotd_25$training)))
                               
for (h1 in unique(plotd_25$host)){
  for (t2 in unique(plotd_25$training)){
#if (h1=="Streptococcus" & t2=="0.25"){next}
wilcx_p[h1,t2] <- try(as.numeric(wilcox.test(x = plotd_25$value[plotd_25$host==h1 &
                                     plotd_25$training==t2 &
                                     plotd_25$stat=="ET"],
                y = plotd_25$value[plotd_25$host==h1 &
                                     plotd_25$training==t2 &
                                     plotd_25$stat=="FPG"])$p.value),
                TRUE)
cat("host",h1,", train",t2,":",wilcx_p[h1,t2],"\n")
  }}

knitr::kable(wilcx_p,digits = 3,format = "latex")

sum(p.adjust(wilcx_p,"holm")<0.05)/(nrow(wilcx_p)*ncol(wilcx_p))
