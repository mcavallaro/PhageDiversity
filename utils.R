#' Function to get species counts from list of individuals observed
#' argument: individuals
#' spec_byhost_l[['Escherichia']]  %>% getSpeciesCount %>%  head()
getSpeciesCount<-function(individuals){
  return(c(table(individuals)))
}

#' Function to get frequency table from  species count
#' argument: speciesCount (table)
getFrequencyTable<-function(speciesCount){
  return(c(table(speciesCount)))
}

#' Function to extract Mr,n from species counts
#' argument: species counts (table)
#M = spec_byhost_l[['Escherichia']] %>% unlist() %>% table%>% extractM()
#M = spec_byhost_l[['Escherichia']] %>% getSpeciesCount %>% extractM()
extractM<-function(counts){
  n<-sum(counts) #No of sampled indiv
  M<-rep(0,n) #at least count 0 in each class
  # for each species, +1 the right class of M
  for (i in seq(along=counts)){
    M[counts[i]]<-M[counts[i]]+1
  }
  return(M)
}

bootstrapSpeciesCount<-function(speciesCount, n_bt_samples=2){
  individuals<-unlist(speciesCount)
  len = length(speciesCount)
  tmp = list()
  for (i in 1:n_bt_samples){
    bootstrapped_speccounts<-sample(speciesCount, len, replace=T)
    tmp[[i]]<-bootstrapped_speccounts
  }
  return(tmp)
}

bootstrapObservedIndividuals<-function(individuals, n_bt_samples=2){
  # this function doesn't make really sense, since it always decreases the diversity
  individuals<-unlist(individuals)
  len<-length(individuals)
  # speccounts <- as.vector(unname(table(individuals)))
  # freq = c(table(speccounts))
  tmp<-list()
  for (i in 1:n_bt_samples){
    bootstrapped_individuals = sample(individuals, len, replace=T)
    # bootstrapped_speccounts = as.vector(unname(table(bootstrapped_individuals)))
    # bootstrapped_freq = c(table(bootstrapped_speccounts))
    tmp[[i]]<-bootstrapped_individuals
  }
  # list(freq_table=freq, bootstrapped_individuals=tmp)
  return(tmp)
}

BootstrapPredictionIntervals<-function(new_species_estimate, m, num_boostrap_samples){
  n = max(m)
  for(i in 1:num_boostrap_samples){
    attribute_name = sprintf("Boostrap sample %d", i)
    n1 = length(attributes(new_species_estimate)[[attribute_name]])
    if (n1 < n){
      n = n1
    }
  }
  tmp.df = data.frame(1:n)
  for(i in 1:num_boostrap_samples){
    attribute_name = sprintf("Boostrap sample %d", i)
    tmp.df$attribute_name = attributes(new_species_estimate)[[attribute_name]][1:n]
    names(tmp.df)[names(tmp.df) == "attribute_name"] <- attribute_name
  }
  lower = apply(tmp.df, 1, function(x){quantile(x, 0.025)})
  upper = apply(tmp.df, 1, function(x){quantile(x, 0.975)})
  return(data.frame(m=1:n, lower=lower, upper=upper))
}