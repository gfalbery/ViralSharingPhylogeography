# Trying different viral sharing models using DNA and RNA viruses ####

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

FinalHostMatrix[,"RNA"] <- RNAHostdf[with(FinalHostMatrix, paste(Sp,Sp2)), "RNA"]
FinalHostMatrix[,"DNA"] <- DNAHostdf[with(FinalHostMatrix, paste(Sp,Sp2)), "DNA"]

RNADNAComp <- gather(FinalHostMatrix, key = "Virus.Type", value = "VirLink", RNA, DNA)
#
#list(
#  ggplot(RNADNAComp, aes(Space, VirLink, colour = Virus.Type)) + 
#    geom_point(alpha = 0.2) + coord_fixed() +
#    geom_smooth(method = lm, fill = NA) + 
#    stat_smooth(method = lm, geom = "ribbon", fill = NA, lty = 2) +
#    theme(legend.position = "top"),
#  
#  ggplot(RNADNAComp, aes(Phylo2, VirLink, colour = Virus.Type)) + 
#    geom_point(alpha = 0.2) + coord_fixed() +
#    geom_smooth(method = lm, fill = NA) +
#    stat_smooth(method = lm, geom = "ribbon", fill = NA, lty = 2) +
#    theme(legend.position = "top")
#) %>% arrange_ggplot2(ncol = 2)

# Trying different viral sharing models using DNA and RNA viruses ####

library(igraph); library(tidyverse)

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

FinalHostMatrix[,"Vector"] <- VectorHostdf[with(FinalHostMatrix, paste(Sp,Sp2)), "Vector"]

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

FinalHostMatrix[,"NVector"] <- NVectorHostdf[with(FinalHostMatrix, paste(Sp,Sp2)), "NVector"]

VectorComp <- gather(FinalHostMatrix, key = "Virus.Type", value = "VirLink", Vector, NVector)


list(
  ggplot(VectorComp, aes(Space, VirLink, colour = Virus.Type)) + 
    geom_point(alpha = 0.2) + coord_fixed() +
    geom_smooth(method = lm, fill = NA) + 
    stat_smooth(method = lm, geom = "ribbon", fill = NA, lty = 2) +
    theme(legend.position = "top"),
  
  ggplot(VectorComp, aes(Phylo2, VirLink, colour =  Virus.Type)) + 
    geom_point(alpha = 0.2) + coord_fixed() +
    geom_smooth(method = lm, fill = NA) +
    stat_smooth(method = lm, geom = "ribbon", fill = NA, lty = 2) +
    theme(legend.position = "top")
  
) %>% arrange_ggplot2(ncol = 2)

# Lumping sharing DNA-borne viruses in #####

df3 <- rbind(RNADNAComp[RNADNAComp$Virus.Type=="DNA",c("Space","Phylo2","Virus.Type","VirLink")],
             VectorComp[,c("Space","Phylo2","Virus.Type","VirLink")])

list(
  ggplot(df3, aes(Space, VirLink, colour = Virus.Type)) + 
    geom_point(alpha = 0.2) + coord_fixed() +
    geom_smooth(method = lm, fill = NA) + 
    stat_smooth(method = lm, geom = "ribbon", fill = NA, lty = 2) +
    theme(legend.position = "top"),
  
  ggplot(df3, aes(Phylo2, VirLink, colour =  Virus.Type)) + 
    geom_point(alpha = 0.2) + coord_fixed() +
    geom_smooth(method = lm, fill = NA) +
    stat_smooth(method = lm, geom = "ribbon", fill = NA, lty = 2) +
    theme(legend.position = "top")
  
) %>% arrange_ggplot2(ncol = 2)






