# Creating viral subsets to elaborate on sharing patterns ####

# DNA and RNA viruses ####

SubResps <- c("RNA", "DNA", "Vector","NVector")

library(igraph); library(tidyverse); library(ggregplot)

RNAViruses <- VirusTraits %>% filter(vDNAoRNA == "RNA") %>% select(vVirusNameCorrected) %>% unlist %>% intersect(rownames(M))
DNAViruses <- VirusTraits %>% filter(vDNAoRNA == "DNA") %>% select(vVirusNameCorrected) %>% unlist %>% intersect(rownames(M))

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
                                                 by = c("Sp","Sp2"))

FinalHostMatrix <- FinalHostMatrix %>% left_join(DNAHostdf,
                                                 by = c("Sp","Sp2"))

# Vector-Borne ####

VectorViruses <- VirusTraits %>% filter(vDNAoRNA == "RNA"&vVectorYNna == "Y") %>% 
  dplyr::select(vVirusNameCorrected) %>% 
  unlist %>% intersect(rownames(M))
MVector <- M[VectorViruses,]
MVector <- MVector[,which(colSums(MVector)>0)]
VectorBipGraph <- graph.incidence(MVector, weighted = T)
VectorHostGraph <- bipartite.projection(VectorBipGraph)$proj2
VectorHostAdj <- get.adjacency(VectorHostGraph) %>% as.matrix

VectorHostdf <- VectorHostAdj %>% reshape2::melt() %>%
  #dplyr::rename(Sp = Var2, Sp2 = Var1, Vector = value) %>%
  dplyr::select(Var2, Var1, value) %>%
  dplyr::rename(Sp = Var2, Sp2 = Var1, Vector = value)

rownames(VectorHostdf) <- with(VectorHostdf, paste(Sp, Sp2))

# Non-Vector-Borne ####

NVectorViruses <- VirusTraits %>% filter(vDNAoRNA == "RNA"&vVectorYNna == "N") %>% select(vVirusNameCorrected) %>% unlist %>% intersect(rownames(M))
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
                                                 by = c("Sp","Sp2"))

FinalHostMatrix <- FinalHostMatrix %>% left_join(NVectorHostdf,
                                                 by = c("Sp","Sp2"))


SlopeTime <- gather(FinalHostMatrix, key = "Group", value = "Shared", paste0(SubResps)) %>%
  filter(!is.na(Shared))

#ggplot(SlopeTime, aes(DietSim, Shared, colour = Group, lty = as.factor(Eaten))) + 
#  facet_wrap(~Group) + 
#  geom_point() + 
#  theme(legend.position = "none") +
#  geom_smooth()

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

