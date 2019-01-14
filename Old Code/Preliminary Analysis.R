
# Preliminary analysis results ####

library(igraph); library(ggregplot); library(ggplot2);library(dplyr)

# Plot virus graph
plot(Virusgraph, vertex.label = NA)

# Virus clustering 

mean_distance(Virusgraph, directed = FALSE)

Viruses[,c("vCytoReplicTF","vSegmentedTF")] <- apply(Viruses[,c("vCytoReplicTF","vSegmentedTF")], 1, as.numeric)

VirCommTests <- c("vGenus","vFamily","vOrder",
                  "vCytoReplicTF",
                  #"vSegmentedTF","vVectorYNna",
                  "vSSoDS","vEnvelope",
                  "vDNAoRNA",
                  "Wildlife","Domestic","Human","HumDomWild"
)

sapply(VirCommTests, function(a){
  igraph::modularity(Virusgraph, as.factor(Viruses[,a]))
})

# Only domestic and non-domestic wildlife viruses are well-clustered
# No phylogenetic clustering??

# Are domestic viruses more central? ####

BarGraph(Viruses, "Domestic", "vDegree", text = "N")
BarGraph(Viruses, "Domestic", "vEigenvector", text = "N")

# Only Eigenvector centrality? 

# Are zoonotic viruses more central?

BarGraph(Viruses, "Human", "vDegree", text = "N")
BarGraph(Viruses, "Human", "vEigenvector", text = "N")

# Nope

BarGraph(Viruses, "HumDomWild", "vDegree", text = "N") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

BarGraph(Viruses, "HumDomWild", "vEigenvector", text = "N") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# human + domestic viruses are most central, but not by much

# Examining clustering by domestic/wildlife ####

Virusgraph2 <- Virusgraph

V(Virusgraph2)$color <- adjustcolor(c("black", "red")[Viruses$Domestic+1], alpha = 0.6)
V(Virusgraph2)$size <- c(3,5)[Viruses$Core+1]

E(Virusgraph2)$color <- adjustcolor("grey", alpha = 0.1)

VirusLayout <- layout_with_kk(Virusgraph2)

plot(Virusgraph2,  # Concentric rings of domestic and non-domestic viruses ####
     layout = VirusLayout, 
     vertex.label = NA)

# Do viruses that infect domestics more often infect humans?

BarGraph(Viruses, "Domestic", "Human")
BarGraph(Viruses[!Viruses$HumanDist==Inf,], 
         "Domestic", "HumanDist",
         text = "N")

# No, domestic viruses *less* frequently infect humans, and further away in the network.
# Are domestic viruses closer to humans in the network, examining all links?

ggMMplot(vHumanDist2, "Domestic", "HumanDist") +
  scale_fill_discrete(limits = c(0:2,"Inf")) + 
  labs(x = "Infects Domestics", 
       y = "vDomestic Distance from Humans",
       title = "Domestic Distances to Humans") + 
  coord_fixed() 

# Domestic viruses are less often connected to humans, 
# but all of them are in some way connected (no infinite distances). 

# Are zoonoses more likely to infect domestic animals?

BarGraph(Viruses, "Human", "Domestic", text = "N")

# Nope

# Is the minimum distance to a domestic animal higher or lower in zoonoses?

BarGraph(Viruses[!NARows(Viruses, "MinDomDist"),], 
         "Human", "MinDomDist", text = "N")

# Human viruses are slightly further away from the nearest domestic animal

# Is the mean distance to a domestic animal higher or lower in human viruses?

BarGraph(Viruses[!NARows(Viruses, "MinDomDist"),], 
         "Human", "MeanDomDist", text = "N")

# Zoonoses' mean distance to domestic animals is much lower

# Do viral families vary in terms of their zoonotic potential?

BarGraph(Viruses[!NARows(Viruses, "HumanDist"),], "vFamily", "HumanDist", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Family", 
       y = "Min Distance to Humans",
       title = "vFamily Distance to Humans") + 
  coord_fixed(ratio = nunique(Viruses$vFamily)) +
  theme(legend.position = "none") +
  ggsave("Figures/vFamily Distance to Humans.tiff", units = "mm", width = 150, height = 150, dpi = 300)

# Yep!

# Do viral families vary in terms of their domestonotic potential?

BarGraph(Viruses[!NARows(Viruses, "MeanDomDist"),], "vFamily", "MeanDomDist", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Family", 
       y = "Min Distance to Domestics",
       title = "vFamily Distance to Domestics") + 
  coord_fixed(ratio = nunique(Viruses$vFamily)/2) +
  theme(legend.position = "none") +
  ggsave("Figures/vFamily Distance to Domestics.tiff", units = "mm", width = 150, height = 150, dpi = 300)

# Yep!

# Does family-level distance from domestics correlate with distance from humans?
BarBarGraph(Viruses[!NARows(Viruses, c("MinDomDist", "HumanDist")),], 
            "vFamily", "MinDomDist", "HumanDist") + 
  theme(legend.position = "none")

# Maybe a negative correlation

BarBarGraph(Viruses[!NARows(Viruses, c("MinDomDist", "HumanDist")),], 
            "vFamily", "MaxDomDist", "HumanDist") + 
  theme(legend.position = "none") + geom_smooth(method = lm)

# But a positive correlation with maximum distance

BarBarGraph(Viruses[!NARows(Viruses, c("MinDomDist", "HumanDist")),], 
            "vFamily", "MeanDomDist", "HumanDist") + 
  theme(legend.position = "none") + geom_smooth(method = lm)

# And with mean distance

# Does species-level distance from domestics correlate with distance from humans?
ggplot(Viruses[!NARows(Viruses, c("MinDomDist", "HumanDist")),], 
       aes(MinDomDist, HumanDist)) + geom_point(position = position_jitter(h = 0.25, w = 0.25)) + 
  geom_smooth()

# Maybe a negative correlation

ggplot(Viruses[!NARows(Viruses, c("MinDomDist", "HumanDist")),], 
       aes(MaxDomDist, HumanDist)) + geom_point(position = position_jitter(h = 0.25, w = 0.25)) + 
  geom_smooth()

# Definitely a positive correlation with with max distance

ggplot(Viruses[!NARows(Viruses, c("MinDomDist", "HumanDist")),], 
       aes(MeanDomDist, HumanDist))  + geom_point(position = position_jitter(h = 0.25, w = 0.25)) + 
  geom_smooth()

# And with mean distance

# To summarise: 
# Domestic viruses are less likely to infect humans
# And they tend to be further away in the network than non-domestic viruses
# BUT they are all in some way connected to humans.

# Minimum distance from domestics correlates negatively with min distance from humans
# BUT mean and maximum distances correlate positively.
# This is true both at the species- and family-level.
# These correlations are probably highly constrained by host breadth etc.
# Also the fact that humans are connected to every single domestic animal,
# so connection to one domestic immediately means a distance of =<2 links to humans

# --> need to start controlling for more biases.


