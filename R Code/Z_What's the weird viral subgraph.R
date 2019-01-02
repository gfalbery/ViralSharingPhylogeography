
# Investigating weird viral subgraph ####

VirusSubgraph1 <- induced_subgraph(Virusgraph, 
                                   Viruses$Degree<200)

plot(VirusSubgraph1, vertex.label = NULL)

VirusSubgraph2 <- induced_subgraph(Virusgraph, 
                                   Viruses$Degree>200)

plot(VirusSubgraph2, vertex.label = NULL)
