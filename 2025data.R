
# import data
fulltable <- read.csv("data/phagesspeciescounts_perhostspec_Sept2024.csv", check.names=F)
spec_byhost <- fulltable |> select(Host, `Phage Species`) |> nest_by(Host)
spec_byhost_l <- as.list(spec_byhost$data)
names(spec_byhost_l) <- spec_byhost$Host
#' count number of samples per host species
#' which ones are > 1000?
nosamples_host <- sapply(spec_byhost_l, nrow)
nospecies_host <- sapply(spec_byhost_l, function(x){nrow(unique(x))})
samples_1K <- which(nosamples_host>999)

# create and print a summary for the hosts with most phage samples
Summary_Sept2024 = data.frame(
  host = c("Escherichia", "Klebsiella", "Mycobacterium", "Pseudomonas", "Salmonella", "Vibrio")
)
Summary_Sept2024$number.of.samples = nosamples_host[Summary_Sept2024$host]
Summary_Sept2024$number.of.species = nospecies_host[Summary_Sept2024$host]
print(Summary_Sept2024)

fulltable2025 <- read.csv("data/3May2025_data.tsv", check.names=F, sep = '\t')
spec_byhost2025 <- fulltable2025 |> select(Host, `vOTU`) |> nest_by(Host)
spec_byhost2025_l <- as.list(spec_byhost2025$data)
names(spec_byhost2025_l) <- spec_byhost2025$Host
nosamples_host2025 <- sapply(spec_byhost2025_l, nrow)
nospecies_host2025 <- sapply(spec_byhost2025_l, function(x){nrow(unique(x))})
Summary_Sept2025 = data.frame(
  host = c("Escherichia", "Klebsiella", "Mycobacterium", "Pseudomonas", "Salmonella", "Vibrio")
)
Summary_Sept2025$number.of.samples = nosamples_host2025[Summary_Sept2025$host]
Summary_Sept2025$number.of.species = nospecies_host2025[Summary_Sept2025$host]
print(Summary_Sept2025)

Summary=merge(Summary_Sept2024, Summary_Sept2025, by='host', suffixes=c(".2024",".2025"))
Summary$new_samples = Summary$number.of.samples.2025 - Summary$number.of.samples.2024
Summary$new_species = Summary$number.of.species.2025 - Summary$number.of.species.2024

print(Summary)


addNewDataPoints<-function(Host, summary.df=Summary, col='black'){
  idx = summary.df$host == Host
  x = summary.df$new_samples[idx]
  y = summary.df$new_species[idx]
  points(x, y, col=col)  
}
