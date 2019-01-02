
# 1b_Comparing network-level stuff ####

library(igraph); library(ggforce); library(ggregplot)

# Viruses ####

# X_Community detection ####

components(Virusgraph)$no

# This takes a while ####
VirusCommunities <- cluster_edge_betweenness(Virusgraph, bridges = T)
VirusBridges <- VirusCommunities$bridges

E(Virusgraph)$color <- "black"
#E(Virusgraph)$color[1:length(E(Virusgraph)) %in% VirusBridges] <- "red"

V(Virusgraph)$color <- membership(VirusCommunities)
#V(Virusgraph)$size <- 5

l <- layout_with_kk(Virusgraph)
plot(Virusgraph, layout = l)
plot(Virusgraph, layout = layout_with_mds)

plot(VirusCommunities, Virusgraph, rescale = F, layout = l*6)
plot(VirusCommunities, Virusgraph, layout = l)

#VirusCliques <- cliques(Virusgraph)

VirusCommunitylist <- list()

for(x in 1:max(VirusCommunities$membership)){
  Communitylist[[x]] <- induced_subgraph(Virusgraph, 
                                         VirusCommunities$membership == x)
}

# Community testing ####

mean_distance(Virusgraph, directed = FALSE)

Viruses[,c("vCytoReplicTF","vSegmentedTF")] <- apply(Viruses[,c("vCytoReplicTF","vSegmentedTF")], 1, as.numeric)

VirCommTests <- c("vGenus","vFamily","vOrder",
                  "vCytoReplicTF",
                  #"vSegmentedTF","vVectorYNna",
                  "vSSoDS","vEnvelope",
                  "vDNAoRNA",
                  "Wildlife","Domestic","Human"
)

sapply(VirCommTests, function(a){
  igraph::modularity(Virusgraph, as.factor(Viruses[,a]))
})

# Modularity for domestic viruses is actually the highest!

Virusgraph2 <- Virusgraph

V(Virusgraph2)$color <- adjustcolor(c("black", "red")[Viruses$Domestic+1], alpha = 0.6)
V(Virusgraph2)$size <- 3

E(Virusgraph2)$color <- adjustcolor("grey", alpha = 0.1)

VirusLayout <- layout_with_kk(Virusgraph2)

plot(Virusgraph2,  # Concentric rings of wildlife and domestic viruses ####
     layout = VirusLayout, 
     vertex.label = NA)

pdf("Virusgraph2.pdf", width = 25, height = 25)
plot(Virusgraph2,  # Concentric rings of wildlife and domestic viruses ####
     layout = VirusLayout, 
     vertex.label = NA)
dev.off()

# Investigating Dom/Wildlife subgraphs ####

DomVirusgraph <- induced_subgraph(Virusgraph, 
                                  Viruses$Domestic == 1)

DomVirusLayout <- layout_with_kk(DomVirusgraph)
plot(DomVirusgraph, layout = DomVirusLayout)

pdf("DomViruses.pdf", width = 25, height = 25)
plot(DomVirusgraph, layout = DomVirusLayout*5)
dev.off()

WildlifeVirusgraph <- induced_subgraph(Virusgraph, 
                                       Viruses$Wildlife == 1)

WildlifeVirusLayout <- layout_with_kk(WildlifeVirusgraph)
plot(WildlifeVirusgraph, layout = WildlifeVirusLayout)

pdf("WildlifeViruses.pdf", width = 25, height = 25)
plot(WildlifeVirusgraph, layout = WildlifeVirusLayout*5)
dev.off()

# Making phylogenetic subgraphs ####

vOrderlist <- vFamilylist <- list()

for(x in 1:max(as.numeric(Viruses$vOrder))){
  vOrderlist[[x]] <- induced_subgraph(Virusgraph, 
                                      as.numeric(Viruses$vOrder) == x)
}

for(x in 1:max(as.numeric(Viruses$vFamily))){
  vFamilylist[[x]] <- induced_subgraph(Virusgraph, 
                                       as.numeric(Viruses$vFamily) == x)
}
