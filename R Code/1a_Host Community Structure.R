
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

V(Hostgraph)$color <- membership(HostCommunities)
V(Hostgraph)$size <- 5
#V(Hostgraph)$label <- NULL

vl <- layout_with_kk(Hostgraph)
plot(Hostgraph, layout = vl)
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

sapply(hCommTests, function(a){
  modularity(Hostgraph, as.factor(Hosts[,a]))
})

hOrderlist <- hFamilylist <- list()

for(x in 1:nunique(Hosts$hOrder)){
  hOrderlist[[x]] <- induced_subgraph(Hostgraph, 
                                     as.numeric(as.factor(Hosts$hOrder)) == x)
}

for(x in 1:nunique(Hosts$hFamily)){
  hFamilylist[[x]] <- induced_subgraph(Hostgraph, 
                                     as.numeric(as.factor(Hosts$hFamily)) == x)
}

# Distances from other types of Hosts ####

# dist.from.NYT # in case trying to google/evernote search the original code

hDomDist <- distances(Hostgraph, v = V(Hostgraph),
                           to = V(Hostgraph)[Hosts$hDom == "domestic"], weights=NA)

hHumanDist <- distances(Hostgraph, v = V(Hostgraph),
                           to = V(Hostgraph)[Hosts$hDom == "human"], weights=NA)

hHumanDist2 <- as.data.frame(hHumanDist)
hHumanDist2$hDom <- Hosts$hDom

hHumanDist2 <- hHumanDist2[!hHumanDist2$Homo_sapiens %in% c("Inf",0),]
BarGraph(hHumanDist2, "hDom", "Homo_sapiens")


