# Data import for EHA ###

remove(list = ls())

library(igraph); library(magrittr); library(dplyr); library(ggplot2); require(RCurl); library(readr)

AssocsBase <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% data.frame()
HostTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/hosts.csv") %>% data.frame()
VirusTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/viruses.csv") %>% data.frame()

names(AssocsBase)[1:2] <- c("Virus", "Host")
AssocsBase <-mutate(AssocsBase,
                    Virus = as.factor(Virus),
                    Host = as.factor(Host))

AssocsBase2 <- AssocsBase

# Rabies and humans are both highly central in the networks
#AssocsBase2 <- droplevels(AssocsBase[!AssocsBase$Host == "Homo_sapiens"&
#                               !AssocsBase$Virus == "Rabies_virus",])
AssocsBase2 <- droplevels(AssocsBase[#!AssocsBase$Host == "Homo_sapiens"&
                               !AssocsBase$Virus == "Rabies_virus",])

# Making bipartite projections ####

AssocsTraits <- AssocsBase2[,1:2]

m <- table(AssocsTraits)
M <- as.matrix(m)

bipgraph <- graph.incidence(M, weighted = T)

Virusgraph <- bipartite.projection(bipgraph)$proj1
Hostgraph <- bipartite.projection(bipgraph)$proj2

Virusadj <- as.matrix(get.adjacency(Virusgraph))
Hostadj <- as.matrix(get.adjacency(Hostgraph))

# Deriving metrics from the networks ####

Hosts <- data.frame(Sp = unique(AssocsBase2$Host),
                    Degree = colSums(Hostadj, na.rm = T),
                    Eigenvector = eigen_centrality(Hostgraph)$vector,
                    Kcore = coreness(Hostgraph),
                    Between = betweenness(Hostgraph))

Viruses <- data.frame(Sp = unique(AssocsBase2$Virus),
                      Degree = colSums(Virusadj, na.rm = T),
                      Eigenvector = eigen_centrality(Virusgraph)$vector,
                      Kcore = coreness(Virusgraph),
                      Between = betweenness(Virusgraph))

#GGally::ggpairs(Hosts[!Hosts$Sp=="Homo_sapiens",2:5], lower = list(continuous = "smooth"))
#GGally::ggpairs(Viruses[,2:5], lower = list(continuous = "smooth"))

Hosts <- merge(Hosts, HostTraits, by.x = "Sp", by.y = "hHostNameFinal", all.x = T)
Viruses <- merge(Viruses, VirusTraits, by.x = "Sp", by.y = "vVirusNameCorrected", all.x = T)

names(Hosts)[which(names(Hosts)=="hWildDomFAO")] <- "hDom"

Domestics <- Hosts[Hosts$hDom == "domestic", "Sp"]
Wildlife <- Hosts[Hosts$hDom == "wild", "Sp"]

DomesticViruses <- as.factor(AssocsBase[AssocsBase$Host %in% Domestics, "Virus"])
WildlifeViruses <- as.factor(AssocsBase[AssocsBase$Host %in% Wildlife, "Virus"])

AssocsTraits <- merge(AssocsTraits, HostTraits, by.x = "Host", by.y = "hHostNameFinal", all.x = T)
AssocsTraits <- merge(AssocsTraits, VirusTraits, by.x = "Virus", by.y = "vVirusNameCorrected", all.x = T)

AssocsTraits$Domestic <- ifelse(AssocsTraits$Host%in%Domestics,1,0)
AssocsTraits$Wildlife <- ifelse(AssocsTraits$Host%in%Wildlife,1,0)

Viruses <- Viruses %>%
  mutate(
    Human = IsZoonotic,
    
    Domestic = case_when(
      Sp %in% DomesticViruses ~ 1,
      TRUE ~ 0),
    
    Wildlife = case_when(
      Sp %in% WildlifeViruses ~ 1,
      TRUE ~ 0),
    
    DomesticCount = c(table(AssocsTraits[AssocsTraits$Domestic == 1,"Virus"])),
    WildlifeCount = c(table(AssocsTraits[AssocsTraits$Wildlife == 1,"Virus"])),
    
    PropDomestic = c(table(AssocsTraits[AssocsTraits$Domestic == 1, "Virus"])/
                       table(AssocsTraits$Virus)),
    
    Records = c(table(AssocsTraits$Virus))
  ) %>%
  select(-IsZoonotic)

vCentrality <- c("vDegree", "vEigenvector", "vCore")
vDists <- c("gaussian", "beta", "binomial")

Viruses$vDegree <- kader:::cuberoot(Viruses$Degree)
Viruses$vEigenvector <- Viruses$Eigenvector
Viruses$vCore <- ifelse(Viruses$Kcore==max(Viruses$Kcore),1,0)

Viruses$vEigenvector[round(Viruses$vEigenvector,4)==0] <- 0.0001 # Needed for a beta-distribution
Viruses$vEigenvector[round(Viruses$vEigenvector,4)==1] <- 0.999

Viruses$vGenomeAveLengthLn <- log(Viruses$vGenomeAveLength)
Viruses$vPubMedCitesLn <- log(Viruses$vPubMedCites + 1)

# Loading functions, determining themes ####

#devtools::install_github("gfalbery/ggregplot")
library(ggregplot)

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

AlberTheme <- theme(axis.text.x = element_text(size = 12, colour = "black"), 
               axis.text.y = element_text(size = 12, colour = "black")) + 
  theme(axis.title.x = element_text(vjust = -0.35), 
        axis.title.y = element_text(vjust=1.2)) + 
  theme_bw()

theme_set(AlberTheme)

