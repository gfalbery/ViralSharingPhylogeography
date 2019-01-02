
# 1b_Comparing network-level stuff ####

library(igraph)

# Viruses ####

# X_Community detection ####

components(Virusgraph)$no

VirusCommunities <- cluster_edge_betweenness(Virusgraph, bridges = T)

VirusBridges <- VirusCommunities$bridges

plot(Virusgraph)

E(Virusgraph)$color <- "black"
#E(Virusgraph)$color[1:length(E(Virusgraph)) %in% VirusBridges] <- "red"

V(Virusgraph)$color <- membership(VirusCommunities)
#V(Virusgraph)$size <- 5

l <- layout_with_kk(Virusgraph)
plot(Virusgraph, layout = l)
plot(Virusgraph, rescale = F, layout = l*2)
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
V(Virusgraph2)$labels <- NULL
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

# Distances from other types of Viruses ####

vDomDist <- distances(Virusgraph, v = V(Virusgraph),
                      to = V(Virusgraph)[Viruses$Domestic == 1], weights=NA)

vHumanDist <- distances(Virusgraph, v = V(Virusgraph),
                        to = V(Virusgraph)[Viruses$Human == 1], weights=NA)

NARows <-function(df, vars){
  apply(df[,vars], 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

vHumanDist2 <- vHumanDist[!NARows(vHumanDist),]

vHumanDist2 <- data.frame(Sp = rownames(vHumanDist2),
                          Mean = apply(vHumanDist2, 1, function(a) mean(a, na.rm = T)),
                          Min = apply(vHumanDist2, 1, function(a) min(a, na.rm = T)),
                          Max = apply(vHumanDist2, 1, function(a) max(a, na.rm = T)),
                          Domestic = Viruses$Domestic[!NARows(vHumanDist)],
                          Family = Viruses$vFamily[!NARows(vHumanDist)])

BarGraph(vHumanDist2, "Domestic", "Mean") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x=  "Infects Domestics", y = "Mean Distance from Humans")

BarGraph(vHumanDist2, "Family", "Mean") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Mean Distance from Humans")

ggMMplot(vHumanDist2, "Domestic", "Min") +
  scale_fill_discrete(limits = 0:3)

ggMMplot(vHumanDist2, "Domestic", "Max") +
  scale_fill_discrete(limits = 0:3)

# Putting in long format ####

vHumanDistLong <- reshape2::melt(vHumanDist)
vHumanDistLong$Domestic <- Viruses$Domestic

ggplot(vHumanDistLong, aes(as.factor(Domestic), value, colour = as.factor(Domestic))) + 
  ggforce::geom_sina() + 
  labs(x = "Infects Domestics", 
       y = "Mean Distance from Humans")

ggMMplot(vHumanDistLong, "Domestic", "value")
