
# 1a_Comparing network-level stuff ####

library(igraph)

# Hosts ####

# X_Community detection ####

components(Hostgraph)$no

HostCommunities <- cluster_edge_betweenness(Hostgraph, bridges = T)

HostBridges <- HostCommunities$bridges

plot(Hostgraph)

E(Hostgraph)$color <- "black"
#E(Hostgraph)$color[1:length(E(Hostgraph)) %in% HostBridges] <- "red"

V(Hostgraph)$color <- "white"
V(Hostgraph)$size <- 5
V(Hostgraph)$name <- NULL

vl <- layout_with_kk(Hostgraph)
plot(Hostgraph, layout = vl, labels = FALSE)

plot(Hostgraph, rescale = F, layout = vl*2)
plot(Hostgraph, layout = layout_with_mds)

plot(HostCommunities, Hostgraph, rescale = F, layout = vl*6)
plot(HostCommunities, Hostgraph, layout = vl)

#HostCliques <- cliques(Hostgraph)

hCommunitylist <- list()

for(x in 1:max(HostCommunities$membership)){
  hCommunitylist[[x]] <- induced_subgraph(Hostgraph, 
                                         HostCommunities$membership == x)
}

# Community testing ####

mean_distance(Hostgraph, directed = FALSE)
modularity(Hostgraph, HostCommunities)

hCommTests <- c("hFamily", 
               "hOrder",
               "hMarOTerr",
               #"Population_trend",
               "hDom")

sapply(hCommTests, function(a){ # Clusters ~substantially by Family and Order
  modularity(Hostgraph, as.factor(Hosts[,a]))
})

hOrderlist <- hFamilylist <- list()

for(x in unique(Hosts$hOrder)){
  hOrderlist[[x]] <- induced_subgraph(Hostgraph, 
                                     Hosts$hOrder == x)
}

sapply(hOrderlist, transitivity)

for(x in unique(Hosts$hFamily)){
  hFamilylist[[x]] <- induced_subgraph(Hostgraph, 
                                     Hosts$hFamily == x)
}

sapply(hFamilylist, transitivity)
