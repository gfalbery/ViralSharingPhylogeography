# Incorporating RNA virus data ####

library(stringr)

BrierleyRNA <- read.csv("data/BrierleyRNA.csv", header = T)

BrierleyRNA$Species <- as.factor(str_replace_all(BrierleyRNA$Species, " ", "_"))

RemoveUnderscores <- substr(BrierleyRNA$Species,nchar(BrierleyRNA$Species),nchar(BrierleyRNA$Species)) == "_"

BrierleyRNA$Species[RemoveUnderscores] <- 
  substr(BrierleyRNA$Species[RemoveUnderscores], 1, nchar(BrierleyRNA$Species[RemoveUnderscores])-1)

RNAgraph <- induced.subgraph(Virusgraph, 
                             Viruses$vDNAoRNA == "RNA")

RNAViruses <- Viruses[Viruses$vDNAoRNA == "RNA",]
RNAViruses <- merge(RNAViruses, BrierleyRNA, by.x = "Sp", by.y = "Species")

