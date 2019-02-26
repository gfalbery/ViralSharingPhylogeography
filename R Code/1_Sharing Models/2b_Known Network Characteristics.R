
# Rscript "R Code/1_Sharing Models/2b_Known Network Characteristics.R" ####

library(igraph); library(tidyverse); library(ggregplot); library(parallel)

source("R Code/0_Data Import/0a_EHA Data Import.R")

load("Output Files/KnownSimGraphs.Rdata")

AllNetworkStats <- function(graph){
  
  list(
    
    degree(graph),
    
    eigen_centrality(graph),
    
    #betweenness(graph),
    
    #closeness(graph),
    
    transitivity(graph),
    
    components(graph)
    
  ) %>% return
}

print("Observed")
ObsNetwork <- AllNetworkStats(Hostgraph)

print("1!")
PredNetwork1 <- mclapply(SimGraphs1, AllNetworkStats, mc.cores = 10)

print("1b!")
PredNetwork1b <- mclapply(SimGraphs1b, AllNetworkStats, mc.cores = 10)

save(ObsNetwork, PredNetwork1, PredNetwork1b, file = "Output Files/KnownNetworkStats.Rdata")

