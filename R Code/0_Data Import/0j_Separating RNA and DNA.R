# Trying different viral sharing models using DNA and RNA viruses ####

library(igraph); library(tidyverse)

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
  rename(Sp = Var2, Sp2 = Var1, RNA = value) %>%
  select(Sp, Sp2, RNA)

DNAHostdf <- DNAHostAdj %>% reshape2::melt() %>%
  rename(Sp = Var2, Sp2 = Var1, DNA = value) %>%
  select(Sp, Sp2, DNA)

rownames(RNAHostdf) <- with(RNAHostdf, paste(Sp, Sp2))
rownames(DNAHostdf) <- with(DNAHostdf, paste(Sp, Sp2))

FinalHostMatrix[,"RNA"] <- RNAHostdf[with(FinalHostMatrix, paste(Sp,Sp2)), "RNA"]
FinalHostMatrix[,"DNA"] <- DNAHostdf[with(FinalHostMatrix, paste(Sp,Sp2)), "DNA"]

#RNADNAComp <- gather(FinalHostMatrix, key = "Virus.Type", value = "VirLink", RNA, DNA)
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



