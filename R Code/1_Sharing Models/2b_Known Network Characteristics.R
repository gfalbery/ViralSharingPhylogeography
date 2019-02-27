
# Rscript "R Code/1_Sharing Models/2b_Known Network Characteristics.R" ####

library(igraph); library(tidyverse); library(ggregplot); library(parallel); library(SpRanger)

source("R Code/0_Data Import/0a_EHA Data Import.R")

load("Output Files/KnownSimGraphs.Rdata")

if(file.exists("Output Files/KnownNetworkStats.Rdata")) load("Output Files/KnownNetworkStats.Rdata") else{
  
  print("Observed")
  ObsNetwork <- AllNetworkStats(Hostgraph)
  
  print("1!")
  PredNetwork1 <- mclapply(SimGraphs1, AllNetworkStats, mc.cores = 10)
  
  print("1b!")
  PredNetwork1b <- mclapply(SimGraphs1b, AllNetworkStats, mc.cores = 10)
  
  save(ObsNetwork, PredNetwork1, PredNetwork1b, file = "Output Files/KnownNetworkStats.Rdata")
  
}

# Comparison of degree stuff ####

PredDegrees1 <- map(PredNetwork1, "Degree") %>% bind_cols()
PredDegrees1b <- map(PredNetwork1b, "Degree") %>% bind_cols()

PredDegrees <- data.frame(PredDegree1 = apply(PredDegrees1, 1, mean),
                          PredDegree1b = apply(PredDegrees1b, 1, mean),
                          Sp = names(PredNetwork1b[[1]]$Degree))

Hosts <- Hosts %>% select(-c("PredDegree1","PredDegree1b")) %>%
                            left_join(PredDegrees, by = "Sp")


