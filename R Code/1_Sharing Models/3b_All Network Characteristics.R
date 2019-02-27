
# 3b_Getting Whole Network Characteristics ####

library(igraph); library(tidyverse); library(ggregplot)

AllNetworkStats <- function(graph){
  
  library(igraph); library(tidyverse); library(ggregplot)
  
  list(
    
    Degree = degree(graph),
    
    Eigen = eigen_centrality(graph),
    
    #Between = betweenness(graph),
    
    #Closeness = closenness(graph),
    
    Trans = transitivity(graph),
    
    Components = components(graph)
    
    #Prev = Prev(graph)
    
  ) %>% return
}

print("Making AllGraphs")

load("Output Files/AllSimGs.Rdata")

print("Deriving Characteristics")

AllPredNetwork <- mclapply(AllSimGs, AllNetworkStats, mc.cores = 10)

save(AllPredNetwork, file = "Output Files/KnownNetworkStats.Rdata")

AllDegrees <- map(AllPredNetwork, "Degree")
AllPredDegrees <- map(AllPredDegrees, mean)

Hosts <- left_join(Hosts, AllPredDegrees, by.x = "Sp", by.y = "Host", all.x = T)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 <- left_join(Panth1, AllPredDegrees, by.x = "Sp", by.y = "Host", all.x = T)

FullPolygons <- left_join(FullPolygons, AllPredDegrees, by.x = "Sp", by.y = "Host", all.x = T)

# Identifying within-order links ####

Panth1 <- Panth1 %>% slice(order(Panth1$Sp)) %>% 
  filter(Sp%in%V(AllSimGs[[1]])$name) %>% droplevels

hOrderList <- mclapply(1:length(AllSimGs), function(i){
  lapply(levels(Panth1$hOrder), function(x) induced_subgraph(AllSimGs[[i]], 
                                                             as.character(Panth1$hOrder) == x))
})

InDegree = lapply(hOrderList, function(a) lapply(a, function(b) degree(b)) %>% unlist)
InNames <- lapply(hOrderList[[1]], function(a) V(a)$name) %>% unlist
InDegDF <- data.frame(InDegree = InDegree,
                      Host = InNames)

Panth1 <- left_join(Panth1, InDegDF, by.x = "Sp", by.y = "Host", all.x = T)
Panth1$OutDegree <- with(Panth1, AllPredDegree-InDegree)

save(Panth1, file = "Output Files/Panth1.Rdata")
