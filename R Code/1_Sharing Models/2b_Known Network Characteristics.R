
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
                          #PredDegree1b = apply(PredDegrees1b, 1, mean),
                          Sp = names(PredNetwork1[[1]]$Degree))

Hosts <- Hosts %>% # select(-c("PredDegree1")) %>% #,"PredDegree1b")) %>%
  left_join(PredDegrees, by = "Sp")

ggplot(Hosts, aes(Degree, PredDegree1)) + 
  geom_point(alpha = 0.5, colour = AlberColours[2]) +
  coord_fixed() +
  labs(x = "Observed Degree", y = "Predicted Degree", title = "With Random Effects") +
  ggsave("SIFigures/Degree.PredDegree1.jpeg", units = "mm", width = 100, height = 100)

ggplot(Hosts, aes(Degree, PredDegree1b)) + 
  geom_point(alpha = 0.5, colour = AlberColours[2]) +
  coord_fixed() +
  labs(x = "Observed Degree", y = "Predicted Degree", title = "No Random Effects") +
  lims(x = c(0, max(Hosts$Degree, na.rm = T)), y = c(0, max(Hosts$Degree, na.rm = T))) +
  ggsave("SIFigures/Degree.PredDegree1b.jpeg", units = "mm", width = 100, height = 100)

GGally::ggpairs(Hosts %>% select(contains("Degree")), lower = list(continuous = "smooth"))


