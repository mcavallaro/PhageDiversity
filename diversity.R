
library(iNEXT)
library(magrittr)
library(dplyr) #splitting into hosts
library(ggplot2)
source("plotiNEXT.R")
#' Read in species count data,
#' extract list with species count per host
#' species
fulltable<-read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost<-fulltable |> select(Host,
                                   `Phage Species`) |> 
  nest_by(Host)
spec_byhost_l<-as.list(spec_byhost$data)
names(spec_byhost_l)<-spec_byhost$Host


fulltable<-read.csv("data/3May2025_data.tsv", check.names = F, sep = "\t")
spec_byhost<-fulltable |> select(Host,`vOTU`) |> nest_by(Host)
spec_byhost_l_2025<-as.list(spec_byhost$data)
names(spec_byhost_l_2025)<-spec_byhost$Host


host_names= c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")

for (n1 in host_names){
  print(c(n1, spec_byhost_l[[n1]]$`Phage Species` %>% unique %>% length(), spec_byhost_l[[n1]]$`Phage Species` %>% length() ))
}


results<-list()
for (n1 in host_names){
  species_aboundance<-spec_byhost_l[[n1]] %>% table() %>% sort(decreasing = T) %>% as.numeric()
  # # n. of species  
  # species_aboundance %>% length()
  # #  of individuals 
  # species_aboundance %>% sum
  results[[n1]]<-species_aboundance
}


results_2025<-list()
for (n1 in host_names){
  species_aboundance<-spec_byhost_l_2025[[n1]] %>% table() %>% sort(decreasing = T) %>% as.numeric()
  # # n. of species  
  # species_aboundance %>% length()
  # #  of individuals 
  # species_aboundance %>% sum
  results_2025[[n1]]<-species_aboundance
}

#  this iNEXT function calculated the so-called Hill numbers (which are generalised entropies)
#  with order parameters q (q=1 corresponds to the Shannon entropy)
A<-iNEXT(results, q = 0)
A_2025<-iNEXT(results_2025, q = 0)
A_asympt<-iNEXT(results, q = 0, endpoint = 10000, knots = 10)
to_rank<-A_asympt$iNextEst$size_based[A_asympt$iNextEst$size_based$m==10000,]
to_rank[order(to_rank$qD),]
# Assemblage     m        Method Order.q       qD   qD.LCL    qD.UCL        SC    SC.LCL SC.UCL
# 50     Salmonella 10000 Extrapolation       0 392.0698 357.4991  426.6405 0.9999994 0.9999376      1
# 60 Staphylococcus 10000 Extrapolation       0 405.0414 349.8811  460.2018 0.9999954 0.9999626      1
# 30  Mycobacterium 10000 Extrapolation       0 613.6232 580.2238  647.0226 0.9999203 0.9980799      1
# 40    Pseudomonas 10000 Extrapolation       0 640.8830 580.4915  701.2745 0.9995859 0.9987314      1
# 70  Streptococcus 10000 Extrapolation       0 732.7307 657.6533  807.8081 0.9999462 0.9998624      1
# 20     Klebsiella 10000 Extrapolation       0 848.6874 756.6410  940.7338 0.9988635 0.9975968      1
# 80         Vibrio 10000 Extrapolation       0 886.0357 766.1854 1005.8860 0.9992217 0.9982875      1
# 10    Escherichia 10000 Extrapolation       0 914.5734 836.4500  992.6968 0.9986043 0.9963481      1

tmp<-A_2025$iNextEst$size_based
Observed_2025<-tmp[tmp$Method=='Observed',]
tmp<-A$iNextEst$size_based
Observed_2024<-tmp[tmp$Method=='Observed',]


 
make_diversity_plot<-function(outfile, A, Observed_2024, Observed_2025, ylab='Species diversity (Hill coefficient)'){
  png(filename = outfile,
      width = 4, height = 3, units = "in", pointsize = 6,
      bg = "white", res = 1000)
  par(mar = c(5, 4, 4, 2) - 0.2)
  
  plot(qD ~ m, data=rbind(Observed_2024, Observed_2025), type="n", xlab="", ylab="")
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)
  plot.iNEXT2(A, show.legend = F, add=T, SITE=host_names, pch=1, show.main	= F)
  title(ylab = ylab, line=2)
  title(xlab = "Number of individuals", line=2.5)
  
  # we can tell that hosts Escherichia and Klesbiella have higher phage diversity than Mycobacterium and Pseudomonas
  # pch<-(16+1:length(host_names))%%25
  col<-gg_color_hue(length(host_names))
  for (i in 1:length(host_names)){
    points(Observed_2025$m[i],  Observed_2025$qD[i], pch=16, cex=2, lw=2, col=col[i] )
  }
  pch<-(16+1:length(host_names))%%25
  col<-gg_color_hue(length(host_names))
  NAs<-rep(NA, length(col))
  legend_names = c(host_names, "Observed (2024)", "Observed (2025)", "Extrapolated (from 2024)")
  legend("topright", legend = legend_names, col = c(col,"darkgray","darkgray","darkgray"),
         lty = c(NAs, NA, NA, 2),
         lwd = 2,
         pch = c(rep(15, length(col)), 1, 16),
         pt.cex = c(rep(2, length(col)), 2, 2, NA),
         bty="n")
  dev.off()
}


make_diversity_plot("diversity_plot_q0__.png", A, Observed_2024, Observed_2025,
                    ylab="Observed number of species (Hill number of order q=0)")

A_1<-iNEXT(results, q = 1)
A_2025_1<-iNEXT(results_2025, q = 1)
A_1_asympt<-iNEXT(results, q = 1, endpoint = 10000, knots = 10)
to_rank<-A_1_asympt$iNextEst$size_based[A_1_asympt$iNextEst$size_based$m==10000,]
to_rank[order(to_rank$qD),]
# Assemblage     m        Method Order.q       qD   qD.LCL   qD.UCL        SC    SC.LCL    SC.UCL
# 60 Staphylococcus 10000 Extrapolation       1 150.8850 134.7388 167.0312 0.9999954 0.9999791 1.0000000
# 50     Salmonella 10000 Extrapolation       1 161.7703 149.2396 174.3009 0.9999994 0.9998117 1.0000000
# 30  Mycobacterium 10000 Extrapolation       1 195.0297 183.3668 206.6927 0.9999203 0.9986040 1.0000000
# 40    Pseudomonas 10000 Extrapolation       1 214.7418 197.4630 232.0205 0.9995859 0.9988085 1.0000000
# 10    Escherichia 10000 Extrapolation       1 300.8581 277.0368 324.6794 0.9986043 0.9969442 1.0000000
# 20     Klebsiella 10000 Extrapolation       1 318.0339 291.7573 344.3105 0.9988635 0.9977879 0.9999391
# 80         Vibrio 10000 Extrapolation       1 392.8223 356.7118 428.9328 0.9992217 0.9986818 0.9997615

tmp<-A_2025_1$iNextEst$size_based
Observed_2025_1<-tmp[tmp$Method=='Observed',]
tmp<-A_1$iNextEst$size_based
Observed_2024_1<-tmp[tmp$Method=='Observed',]

make_diversity_plot("diversity_plot_q1__.png", A_1, Observed_2024_1, Observed_2025_1,
                    ylab="Hill-Shannon diversity (q=1)")

A_2<-iNEXT(results, q = 2)
A_2025_2 = iNEXT(results_2025, q = 2)
A_2_asympt = iNEXT(results, q = 2, endpoint = 10000, knots = 10)
to_rank<-A_2_asympt$iNextEst$size_based[A_2_asympt$iNextEst$size_based$m==10000,]
to_rank[order(to_rank$qD), ]
# Assemblage     m        Method Order.q        qD    qD.LCL    qD.UCL        SC    SC.LCL    SC.UCL
# 60 Staphylococcus 10000 Extrapolation       2  36.85243  29.12467  44.58019 0.9999954 0.9999589 1.0000000
# 30  Mycobacterium 10000 Extrapolation       2  51.50771  46.96609  56.04933 0.9999203 0.9984019 1.0000000
# 50     Salmonella 10000 Extrapolation       2  53.88432  45.73562  62.03302 0.9999994 0.9999367 1.0000000
# 40    Pseudomonas 10000 Extrapolation       2  56.79410  46.87573  66.71246 0.9995859 0.9988483 1.0000000
# 10    Escherichia 10000 Extrapolation       2  69.05187  60.40140  77.70233 0.9986043 0.9959824 1.0000000
# 20     Klebsiella 10000 Extrapolation       2  76.51061  64.90712  88.11411 0.9988635 0.9977467 0.9999803
# 80         Vibrio 10000 Extrapolation       2  95.51050  73.05721 117.96378 0.9992217 0.9982337 1.0000000
# 70  Streptococcus 10000 Extrapolation       2 239.25367 208.06257 270.44476 0.9999462 0.9998455 1.0000000


tmp<-A_2025_2$iNextEst$size_based
Observed_2025_2<-tmp[tmp$Method=='Observed',]
tmp<-A_2$iNextEst$size_based
Observed_2024_2<-tmp[tmp$Method=='Observed',]

make_diversity_plot("diversity_plot_q2__.png", A_2, Observed_2024_2, Observed_2025_2,
                    ylab="Hill-Simpson diversity (q=2)")


################

name<-host_names[1]
A_2025_<-iNEXT(results_2025[name], q = 0)
tmp_<-A_2025_$iNextEst$size_based
Observed_2025<-tmp_[(tmp_$Method=='Observed')& (tmp_$Assemblage==name),
                    c("Assemblage","m", "qD", "qD.LCL", "qD.UCL")]
A_<-iNEXT(results[name], q = 0, endpoint = Observed_2025$m)
tmp<-A_$iNextEst$size_based
extrapolation<-tmp[(tmp$Method=='Extrapolation') & (tmp$Assemblage==name) & (tmp$m == Observed_2025$m),]
Observed_2024<-tmp[(tmp$Method=='Observed') & (tmp$Assemblage==name), c("m", "qD", "qD.LCL", "qD.UCL")]
extrapolated_v_observed_q0<-cbind(Observed_2024, extrapolation, Observed_2025)
for (name in host_names[-1]){
  A_2025_<-iNEXT(results_2025[name], q = 0)
  tmp_ = A_2025_$iNextEst$size_based
  Observed_2025 = tmp_[(tmp_$Method=='Observed')& (tmp_$Assemblage==name),
                       c("Assemblage","m", "qD", "qD.LCL", "qD.UCL")]
  A_<-iNEXT(results[name], q = 0, endpoint = Observed_2025$m)
  tmp<-A_$iNextEst$size_based
  extrapolation<-tmp[(tmp$Method=='Extrapolation') & (tmp$Assemblage==name) & (tmp$m == Observed_2025$m),]
  Observed_2024<-tmp[(tmp$Method=='Observed') & (tmp$Assemblage==name), c("m", "qD", "qD.LCL", "qD.UCL")]
  extrapolated_v_observed_q0<-rbind(extrapolated_v_observed_q0, cbind(Observed_2024, extrapolation, Observed_2025))
}

name<-host_names[1]
A_2025_<-iNEXT(results_2025[name], q = 1)
tmp_<-A_2025_$iNextEst$size_based
Observed_2025<-tmp_[(tmp_$Method=='Observed')& (tmp_$Assemblage==name),
                     c("Assemblage","m", "qD", "qD.LCL", "qD.UCL")]
A_<-iNEXT(results[name], q = 1, endpoint = Observed_2025$m)
tmp<-A_$iNextEst$size_based
extrapolation<-tmp[(tmp$Method=='Extrapolation') & (tmp$Assemblage==name) & (tmp$m == Observed_2025$m),]
Observed_2024<-tmp[(tmp$Method=='Observed') & (tmp$Assemblage==name), c("m", "qD", "qD.LCL", "qD.UCL")]
extrapolated_v_observed_q1<-cbind(Observed_2024, extrapolation, Observed_2025)
for (name in host_names[-1]){
  A_2025_<-iNEXT(results_2025[name], q = 1)
  tmp_ = A_2025_$iNextEst$size_based
  Observed_2025 = tmp_[(tmp_$Method=='Observed')& (tmp_$Assemblage==name),
                       c("Assemblage","m", "qD", "qD.LCL", "qD.UCL")]
  A_<-iNEXT(results[name], q = 1, endpoint = Observed_2025$m)
  tmp<-A_$iNextEst$size_based
  extrapolation<-tmp[(tmp$Method=='Extrapolation') & (tmp$Assemblage==name) & (tmp$m == Observed_2025$m),]
  Observed_2024<-tmp[(tmp$Method=='Observed') & (tmp$Assemblage==name), c("m", "qD", "qD.LCL", "qD.UCL")]
  extrapolated_v_observed_q1<-rbind(extrapolated_v_observed_q1, cbind(Observed_2024, extrapolation, Observed_2025))
}



name<-host_names[1]
A_2025_<-iNEXT(results_2025[name], q = 2)
tmp_<-A_2025_$iNextEst$size_based
Observed_2025<-tmp_[(tmp_$Method=='Observed')& (tmp_$Assemblage==name),
                     c("Assemblage","m", "qD", "qD.LCL", "qD.UCL")]
A_<-iNEXT(results[name], q = 2, endpoint = Observed_2025$m)
tmp<-A_$iNextEst$size_based
extrapolation<-tmp[(tmp$Method=='Extrapolation') & (tmp$Assemblage==name) & (tmp$m == Observed_2025$m),]
Observed_2024<-tmp[(tmp$Method=='Observed') & (tmp$Assemblage==name), c("m", "qD", "qD.LCL", "qD.UCL")]
extrapolated_v_observed_q2<-cbind(Observed_2024, extrapolation, Observed_2025)
for (name in host_names[-1]){
  A_2025_<-iNEXT(results_2025[name], q=2)
  tmp_<-A_2025_$iNextEst$size_based
  Observed_2025<-tmp_[(tmp_$Method=='Observed')& (tmp_$Assemblage==name),
                       c("Assemblage","m", "qD", "qD.LCL", "qD.UCL")]
  A_<-iNEXT(results[name], q = 2, endpoint = Observed_2025$m)
  tmp<-A_$iNextEst$size_based
  extrapolation<-tmp[(tmp$Method=='Extrapolation') & (tmp$Assemblage==name) & (tmp$m == Observed_2025$m),]
  Observed_2024<-tmp[(tmp$Method=='Observed') & (tmp$Assemblage==name), c("m", "qD", "qD.LCL", "qD.UCL")]
  extrapolated_v_observed_q2<-rbind(extrapolated_v_observed_q2, cbind(Observed_2024, extrapolation, Observed_2025))
}


write.table(extrapolated_v_observed_q0, file = "extrapolated_v_observed_q0.csv", sep = ',')
write.table(extrapolated_v_observed_q1, file = "extrapolated_v_observed_q1.csv", sep = ',')
write.table(extrapolated_v_observed_q2, file = "extrapolated_v_observed_q2.csv", sep = ',')


kable(extrapolated_v_observed_q0[,c(1,2,3,4,6,9,10,11,17,18,19)]%>% round(2), format = 'latex') 
# |    |Assemblage     |    m|       qD|   qD.LCL|   qD.UCL| qD.1| qD.LCL.1| qD.UCL.1|
#   |:---|:--------------|----:|--------:|--------:|--------:|----:|--------:|--------:|
#   |40  |Escherichia    | 2863| 749.4417| 726.1748| 772.7085|  826| 801.8322| 850.1678|
#   |401 |Klebsiella     | 1821| 587.8128| 560.2491| 615.3764|  616| 593.8549| 638.1451|
#   |402 |Mycobacterium  | 3195| 566.2768| 550.3731| 582.1805|  574| 556.5817| 591.4183|
#   |403 |Pseudomonas    | 1878| 476.9020| 456.0594| 497.7446|  517| 491.7363| 542.2637|
#   |404 |Salmonella     | 1584| 346.5091| 331.3074| 361.7108|  370| 356.2474| 383.7526|
#   |405 |Staphylococcus | 1005| 291.6682| 273.2157| 310.1206|  313| 295.9553| 330.0447|
#   |406 |Streptococcus  | 1331| 536.5540| 504.0873| 569.0207|  705| 678.9928| 731.0072|
#   |407 |Vibrio         | 2381| 711.0964| 657.5078| 764.6851|  602| 576.8234| 627.1766|


kable(extrapolated_v_observed_q1[,c(1,2,3,4,6,9,10,11,17,18,19)] %>% round(2), format = 'latex') 
# |    |Assemblage     |    m|       qD|   qD.LCL|   qD.UCL|     qD.1| qD.LCL.1| qD.UCL.1|
#   |:---|:--------------|----:|--------:|--------:|--------:|--------:|--------:|--------:|
#   |40  |Escherichia    | 2863| 253.6587| 237.3079| 270.0095| 268.0111| 252.1937| 283.8284|
#   |401 |Klebsiella     | 1821| 245.0353| 223.9556| 266.1150| 233.6961| 215.5052| 251.8869|
#   |402 |Mycobacterium  | 3195| 176.9624| 168.5902| 185.3346| 176.2516| 166.2090| 186.2942|
#   |403 |Pseudomonas    | 1878| 177.1305| 163.5033| 190.7577| 185.3324| 169.2617| 201.4031|
#   |404 |Salmonella     | 1584| 140.8904| 131.9623| 149.8186| 152.2458| 141.2250| 163.2665|
#   |405 |Staphylococcus | 1005| 120.0394| 107.6366| 132.4422| 125.7000| 114.5058| 136.8942|
#   |406 |Streptococcus  | 1331| 363.0400| 340.5837| 385.4962| 508.6822| 475.6723| 541.6920|
#   |407 |Vibrio         | 2381| 328.5943| 297.8885| 359.3001| 140.5453| 128.3832| 152.7075|
#   
#   
kable(extrapolated_v_observed_q2[,c(1,2,3,4,6,9,10,11,17,18,19)] %>% round(2), format = 'latex') 
# |    |Assemblage     |    m|        qD|    qD.LCL|    qD.UCL|      qD.1|  qD.LCL.1|  qD.UCL.1|
#   |:---|:--------------|----:|---------:|---------:|---------:|---------:|---------:|---------:|
#   |40  |Escherichia    | 2863|  67.89988|  61.05221|  74.74756|  66.49174|  59.00341|  73.98007|
#   |401 |Klebsiella     | 1821|  74.00060|  64.46427|  83.53692|  61.66510|  53.46444|  69.86577|
#   |402 |Mycobacterium  | 3195|  50.95945|  45.38900|  56.52990|  50.15169|  45.54591|  54.75747|
#   |403 |Pseudomonas    | 1878|  55.45582|  46.84309|  64.06855|  55.82809|  48.80781|  62.84838|
#   |404 |Salmonella     | 1584|  52.41151|  46.17790|  58.64511|  58.12575|  48.71851|  67.53298|
#   |405 |Staphylococcus | 1005|  35.70653|  29.18890|  42.22417|  35.32421|  28.54567|  42.10274|
#   |406 |Streptococcus  | 1331| 207.11131| 183.10676| 231.11585| 312.72039| 282.28263| 343.15815|
#   |407 |Vibrio         | 2381|  92.70653|  74.35800| 111.05506|  34.21237|  30.49480|  37.92994|



hills_number<-function(a, q){
  n<-sum(a)
  if(q==1){
    H<-sum(a/n * log(a/n))
    return(exp(-H))
  }else{
    ret<-sum((a / n)**q)
    return(ret ** (1/(1-q)))        
  }
}

