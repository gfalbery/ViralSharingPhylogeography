# Creating viral subsets to elaborate on sharing patterns ####

# DNA and RNA viruses ####

library(igraph); library(tidyverse); library(ggregplot)

RNAViruses <- Viruses %>% filter(vDNAoRNA == "RNA") %>% select(Sp) %>% unlist
DNAViruses <- Viruses %>% filter(vDNAoRNA == "DNA") %>% select(Sp) %>% unlist

MRNA <- M[RNAViruses,]
MDNA <- M[DNAViruses,]

MRNA <- MRNA[,which(colSums(MRNA)>0)]
MDNA <- MDNA[,which(colSums(MDNA)>0)]

RNABipGraph <- graph.incidence(MRNA, weighted = T)
DNABipGraph <- graph.incidence(MDNA, weighted = T)

RNAHostGraph <- bipartite.projection(RNABipGraph)$proj2
DNAHostGraph <- bipartite.projection(DNABipGraph)$proj2

RNAHostAdj <- get.adjacency(RNAHostGraph) %>% as.matrix
DNAHostAdj <- get.adjacency(DNAHostGraph) %>% as.matrix

RNAHostdf <- RNAHostAdj %>% reshape2::melt() %>%
  select(Var2, Var1, value) %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, RNA = value)

DNAHostdf <- DNAHostAdj %>% reshape2::melt() %>%
  select(Var2, Var1, value) %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, DNA = value)

rownames(RNAHostdf) <- with(RNAHostdf, paste(Sp, Sp2))
rownames(DNAHostdf) <- with(DNAHostdf, paste(Sp, Sp2))

FinalHostMatrix <- FinalHostMatrix %>% left_join(RNAHostdf,
                                                 by = c("Sp","Sp2"), all.x = T)

FinalHostMatrix <- FinalHostMatrix %>% left_join(DNAHostdf,
                                                 by = c("Sp","Sp2"), all.x = T)

# Vector-Borne ####

VectorViruses <- Viruses %>% filter(vDNAoRNA == "RNA"&vVectorYNna == "Y") %>% select(Sp) %>% unlist
MVector <- M[VectorViruses,]
MVector <- MVector[,which(colSums(MVector)>0)]
VectorBipGraph <- graph.incidence(MVector, weighted = T)
VectorHostGraph <- bipartite.projection(VectorBipGraph)$proj2
VectorHostAdj <- get.adjacency(VectorHostGraph) %>% as.matrix

VectorHostdf <- VectorHostAdj %>% reshape2::melt() %>%
  #dplyr::rename(Sp = Var2, Sp2 = Var1, Vector = value) %>%
  select(Var2, Var1, value) %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, Vector = value)

rownames(VectorHostdf) <- with(VectorHostdf, paste(Sp, Sp2))

# Non-Vector-Borne ####

NVectorViruses <- Viruses %>% filter(vDNAoRNA == "RNA"&vVectorYNna == "N") %>% select(Sp) %>% unlist
MNVector <- M[NVectorViruses,]
MNVector <- MNVector[,which(colSums(MNVector)>0)]
NVectorBipGraph <- graph.incidence(MNVector, weighted = T)
NVectorHostGraph <- bipartite.projection(NVectorBipGraph)$proj2
NVectorHostAdj <- get.adjacency(NVectorHostGraph) %>% as.matrix

NVectorHostdf <- NVectorHostAdj %>% reshape2::melt() %>%
  select(Var2, Var1, value) %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, NVector = value)

rownames(NVectorHostdf) <- with(NVectorHostdf, paste(Sp, Sp2))

FinalHostMatrix <- FinalHostMatrix %>% left_join(VectorHostdf,
                                                 by = c("Sp","Sp2"), all.x = T)

FinalHostMatrix <- FinalHostMatrix %>% left_join(NVectorHostdf,
                                                 by = c("Sp","Sp2"), all.x = T)



