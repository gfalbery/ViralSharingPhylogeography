
# 3b_Getting Whole Network Characteristics ####

library(igraph); library(tidyverse); library(ggregplot); library(parallel)

source("R Code/00_Master Code.R")

print("Making AllGraphs")

load("Output Files/AllSimGs.Rdata")
load("Output Files/AllSims.Rdata")

print("Deriving Characteristics")

if(file.exists("Output Files/AllNetworkStats.Rdata")) load("Output Files/AllNetworkStats.Rdata") else{
  
  AllPredNetwork <- mclapply(AllSimGs, AllNetworkStats, mc.cores = 10)
  
  save(AllPredNetwork, file = "Output Files/AllNetworkStats.Rdata")
  
}

print("Getting links")

AllDegrees <- map(AllPredNetwork, "Degree") %>% bind_cols()

AllPredDegrees <- data.frame(AllPredDegree = apply(AllDegrees, 1, mean),
                             Sp = rownames(AllSims[[1]]))

Hosts <- left_join(Hosts, AllPredDegrees, by = "Sp", all.x = T)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 <- left_join(Panth1, AllPredDegrees, by.x = "Sp", by.y = "Host", all.x = T)

load("data/FullPolygons.Rdata")

FullPolygons <- left_join(FullPolygons, AllPredDegrees, by = c("Host" = "Sp"), all.x = T)

print("Dividing within- and between-order links")

# Identifying within-order links ####

Panth1 <- Panth1 %>% slice(order(Panth1$Sp)) %>% 
  filter(Sp%in%V(AllSimGs[[1]])$name) %>% droplevels

hOrderList <- mclapply(1:length(AllSimGs), function(i){
  lapply(levels(Panth1$hOrder), function(x) induced_subgraph(AllSimGs[[i]], 
                                                             as.character(Panth1$hOrder) == x))
}, mc.cores = 10)

InDegree <- lapply(hOrderList, function(a) lapply(a, function(b) degree(b)) %>% unlist) %>% bind_cols
InDegrees <- apply(InDegree, 1, mean)
InNames <- lapply(hOrderList[[1]], function(a) V(a)$name) %>% unlist
InDegDF <- data.frame(InDegree = InDegrees,
                      Sp = InNames)

Panth1 <- left_join(Panth1, InDegDF, by = "Sp", all.x = T)
Panth1$OutDegree <- with(Panth1, AllPredDegree-InDegree)

save(Panth1, file = "Output Files/Panth1.Rdata")
