library(knitr)
library(tidyverse)
source("utils.R")
source("import_data.R")

host_names<-c("Escherichia", "Klebsiella",
              "Mycobacterium", "Pseudomonas",
              "Salmonella", "Staphylococcus",
              "Streptococcus", "Vibrio")


fulltable<-read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost<-fulltable |> select(Host,
                                   `Phage Species`) |> 
  nest_by(Host)
spec_byhost_l<-as.list(spec_byhost$data)
names(spec_byhost_l)<-spec_byhost$Host
#' count number of samples per host species
#' which ones are > 1000?
nosamples_host<-sapply(spec_byhost_l, nrow)
samples_1K<-which(nosamples_host > 999)

my_hist<-function(freq_table, titl, ...){
  x<-as.numeric(names(freq_table) )
  y<-freq_table
  X<-c(0.5)
  Y<-c(0)
  for (i in 1:length(y)){
    X<-c(X, x[i]-0.5, x[i]-0.5, x[i]+0.5, x[i]+0.5)
    Y<-c(Y, 0, y[i], y[i], 0)
  }
  plot(X, Y,
	   type = 'n',
	   xaxt = "n",
	   xlab = '',
	   ylab = '',
	   axes = F, ...)
  axis(1, at = x, las = 1, lwd = 0.2)  
  axis(2, las = 1, lwd = 0.2)
  title(ylab = "Frequency", line = 2.6)
  title(xlab = "Phage species count", line = 2, cex = 1.3)  
  title(titl, adj = 0, cex = 2)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.4)
  polygon(X, Y, col = 'lightblue', border = "darkblue", lwd = 0.1)
  usr<-par("usr")  # returns c(x1, x2, y1, y2)
  
  # Calculate width and height of the plot
  plot_width<-usr[2] - usr[1]
  plot_height<-usr[4] - usr[3]
  xa<-1.2
  ya<-1.4
  # Define bottom-left corner of the square
  x0<-usr[2] - plot_width / 2 * xa
  y0<-usr[4] - plot_height / 2 * ya
   
  # Draw the square
  rect(x0, y0, x0 + plot_width / 2  * xa, y0 + plot_height / 2 * ya,
	   border = "white",
	   col="white")
}

make_plot<-function(spec_byhost_l, outfile){
  png(filename = outfile,
      width = 4, height = 4, units = "in", pointsize = 5,
      bg = "white", res = 1000)
  par(mfrow = c(4, 2), mar=c(5, 4, 4, 2) - 0.2)
  for (n1 in host_names){
    print(n1)
    speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
    freq_table<-getFrequencyTable(speccounts)
    my_hist(freq_table, titl=n1)
  }
  
  # Fig = list(
  #   c(0.25, 0.5, 0.68, 1),
  #   c(0.75, 1, 0.68, 1),
  #   c(0.25, 0.5, 0.18, 0.5),
  #   c(0.75, 1, 0.18, 0.5)
  # )
  
  # Fig = list(
  #   c(0.25, 0.5, 0.777, 1),
  #   c(0.75, 1, 0.7777, 1),
  #   c(0.25, 0.5, 0.666666 - 0.2222, 0.6666),
  #   c(0.75, 1, 0.666666 - 0.22222, 0.6666),
  #   c(0.25, 0.5, 0.11111, 0.33333),
  #   c(0.75, 1, 0.11111111, 0.3333)
  # )
  
  Fig<-list(
    c(0.25, 0.5, 1 - 0.1668168, 1),
    c(0.75, 1,   1 - 0.1668168, 1),
    c(0.25, 0.5, 0.75 - 0.1668168, 0.75),
    c(0.75, 1,   0.75 - 0.1668168, 0.75),
    c(0.25, 0.5, 0.5 - 0.1668168, 0.50),
    c(0.75, 1,   0.5 - 0.1668168, 0.50),
    c(0.25, 0.5, 0.25 - 0.1668168, 0.25),
    c(0.75, 1,   0.25 - 0.1668168, 0.25)
  )
  
  for (i in 1:length(host_names)){
    n1<-host_names[i]
    speccounts<-getSpeciesCount(spec_byhost_l[[n1]])
    freq_table<-getFrequencyTable(speccounts)
    x<-as.numeric(names(freq_table) )
    y<-freq_table
    par(fig = unlist(Fig[[i]]), new = TRUE)
    plot(x, y,
	     type = 'n',
		 xaxt = "n",
		 xlab = '',
		 ylab = '',
		 axes = F,
		 log = 'xy',
		 lwd = 0.5)
    axis(1, las = 1, lwd = 0.1)
    axis(2, las = 1, lwd = 0.1)
    lines(x, y, col = 'lightblue', lty = 4, lwd = 0.5)
    points(x, y, col = 'darkblue', pch = 1, cex = 0.7, lwd = 0.5)
    title(ylab = "Frequency", cex = 0.7)
    title(xlab = "Phage species count", line = 2, cex = 0.08)
  }
  dev.off()
}


make_plot(spec_byhost_l, 'hist_2024_.png')


fulltable<-read.csv("data/3May2025_data.tsv", check.names = F, sep = "\t")
spec_byhost<-fulltable |> select(Host,`vOTU`) |> nest_by(Host)
spec_byhost_l_2025<-as.list(spec_byhost$data)
names(spec_byhost_l_2025)<-spec_byhost$Host

make_plot(spec_byhost_l_2025, 'hist_2025_.png')
