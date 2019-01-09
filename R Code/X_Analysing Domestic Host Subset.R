
# Analysing Domestic subset ####

library(igraph); library(ggnet); library(dplyr); library(ggregplot)

DomHostGraph <- induced_subgraph(Hostgraph, Hosts$hDom == "domestic")

plot(DomHostGraph)

DomLayout <- layout_with_kk(DomHostGraph)

ggnet2(DomHostGraph, mode = DomLayout, label = T)

DomFactors <- c("DOMYearBP", "domestic_category", "hContinents")
PhyloFactors <- c("hOrder", "hFamily", "hSpecies")

DomHosts <- Hosts %>%
  filter(hDom == "domestic")

DomHosts[is.na(DomHosts$domestic_category),"domestic_category"] <- "Peridomesticated"

sapply(c(DomFactors, PhyloFactors), function(a){
  igraph::modularity(DomHostGraph, as.factor(DomHosts[,a]))
})

DomHosts %>% mutate(
  dEigenvector <- eigen_centrality(DomHostGraph)$vector,
  dDegree <- rowSums(as.matrix(get.adjacency(DomHostGraph))),
  dKcore <- coreness(DomHostGraph)
)

DomHosts$dCore <- ifelse(DomHosts$Kcore==max(DomHosts$Kcore), 1, 0)

BarGraph(DomHosts, "domestic_category", "dDegree", text = "N") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

BarGraph(DomHosts, "domestic_category", "dEigenvector", text = "N") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




