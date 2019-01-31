# Incorporating RNA virus data ####

library(stringr); library(grid); library(tidyverse); library(igraph); library(ggnet)

# Virus name disputes ####

# Downloading Brierley et al. (2018) data: https://datashare.is.ed.ac.uk/handle/10283/2970
BrierleyRNA <- read.csv("data/BrierleyRNA.csv", header = T)

BrierleyRNA$Species <- str_replace_all(BrierleyRNA$Species, " ", "_")

RemoveUnderscores <- substr(BrierleyRNA$Species,nchar(BrierleyRNA$Species),nchar(BrierleyRNA$Species)) == "_"

BrierleyRNA$Species[RemoveUnderscores] <- 
  substr(BrierleyRNA$Species[RemoveUnderscores], 1, nchar(BrierleyRNA$Species[RemoveUnderscores])-1)

BrierleyRNA <- BrierleyRNA[order(BrierleyRNA$Species),]

RNAViruses <- Viruses[Viruses$vDNAoRNA == "RNA",]

write.csv(setdiff(RNAViruses$Sp, BrierleyRNA$Species), file = "data/RNANames.csv")
write.csv(setdiff(BrierleyRNA$Species, RNAViruses$Sp), file = "data/BrierleyNames.csv")

# By hand, go through to identify synonyms then import set of synonymous pairs of names
# Could probably do something with ICTV lists 
# Liam doesn't have any datasets/lists of synonyms that could help

RNASynonyms <- read.csv("data/RNASynonyms.csv", header = T)

BrierleyRNA$Sp <- sapply(BrierleyRNA$Species, function(a) ifelse(as.character(a)%in%as.character(RNASynonyms$BrierleyName), 
                                                                 as.character(RNASynonyms[as.character(RNASynonyms$BrierleyName)==a, "HP3Name"]),
                                                                 as.character(a)))

# Making the datasets ####

RNAViruses <- Viruses[Viruses$vDNAoRNA == "RNA",]
RNAViruses <- merge(RNAViruses, BrierleyRNA, by.x = "Sp", by.y = "Sp")
RNAViruses <- RNAViruses %>% rename("Iatrogenic" = "Iatrogenic..inc..blood.")

RNAcovar <- c("Vector", "Inhalation", "Ingestion", "Sexual", "Iatrogenic",
              "Fomites", "Broken.Skin", "Maternal", "Direct.Contact",
              "Person.to.person", "Host.range", "Human.only")

for(x in RNAcovar[-which(RNAcovar=="Host.range")]){
  RNAViruses[RNAViruses[,x]=="1*",x] <- 1
  RNAViruses[RNAViruses[,x]=="?",x] <- 0
  RNAViruses[,x] <- as.numeric(RNAViruses[,x])-min(as.numeric(RNAViruses[,x]))
}

RNAViruses <- RNAViruses %>% 
  mutate(
    Host.range = ifelse(RNAViruses$Host.range == "unknown", "narrow", as.character(RNAViruses$Host.range)),
    nTransmission.level = as.numeric(as.factor(Transmission.level)))

FullRNAgraph <- induced.subgraph(Virusgraph, 
                                 Viruses$vDNAoRNA == "RNA")

RNAgraph <- induced.subgraph(Virusgraph, 
                             V(Virusgraph)$name %in% as.character(RNAViruses$Sp))

ggnet2(RNAgraph, label = T)

sapply(RNAcovar, function(a){ # Calculating modularity
  igraph::modularity(RNAgraph, as.factor(RNAViruses[,a]))
}) %>% round(digits = 3)

RNAViruses <- RNAViruses %>% mutate(
  rEigenvector = eigen_centrality(RNAgraph)$vector,
  rDegree = rowSums(as.matrix(get.adjacency(RNAgraph))),
  rKcore = coreness(RNAgraph)
)

# no clustering 

# How do transmission modes etc change with zoonotic/domestonotic status
RNAcovar %>% lapply(function(a){ 
  ggMMplot(RNAViruses, a, "Human") +
    theme(legend.position = "top") +
    ggtitle(a) + coord_fixed()
})%>% arrange_ggplot2(ncol = 4)

RNAcovar %>% lapply(function(a){ 
  ggMMplot(RNAViruses, a, "Domestic") +
    theme(legend.position = "top") +
    ggtitle(a) + coord_fixed()
})%>% arrange_ggplot2(ncol = 4)

# parsing discrepancies
RNAViruses %>% filter(Person.to.person == 1 & Human == 0) # Error or not known in 2017?
RNAViruses %>% filter(Human.only == 1 & Domestic == 1) # Incomplete Brierley data?

# Investigating tripartite associations 
RNAcovar %>% lapply(function(a){ 
  ggMMplot(RNAViruses, a, "HumDomWild") +
    theme(legend.position = "right") +
    ggtitle(a) + coord_fixed()
})%>% arrange_ggplot2(ncol = 4)

ggMMplot(RNAViruses, "Vector", "HumDomWild") + # Most interesting one?
  theme(legend.position = "right") +
  ggtitle("Vector-Borne") + coord_fixed()

RNAViruses %>% filter(Sexual == 1 & Wildlife == 1) %>% select(Sp, Human) # Incomplete Brierley data?

# Looking at transmission level 
ggMMplot(RNAViruses, "Transmission.level", "HumDomWild") + # Most interesting one?
  theme(legend.position = "right") +
  coord_fixed()

RNAcovar %>% lapply(function(a){ 
  BarGraph(RNAViruses, a, "nTransmission.level", text = "N") +
    theme(legend.position = "none") +
    ggtitle(a) + coord_fixed()
})%>% arrange_ggplot2(ncol = 4)

# Looking at centrality
RNAcovar %>% lapply(function(a){ 
  BarGraph(RNAViruses, a, "rDegree", text = "N") +
    theme(legend.position = "none") +
    ggtitle(a) #+ coord_fixed()
})%>% arrange_ggplot2(ncol = 4)

RNAcovar %>% lapply(function(a){ 
  BarGraph(RNAViruses, a, "rEigenvector", text = "N") +
    theme(legend.position = "none") +
    ggtitle(a) #+ coord_fixed()
})%>% arrange_ggplot2(ncol = 4)

# Looking at 


