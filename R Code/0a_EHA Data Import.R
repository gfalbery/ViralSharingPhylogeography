# Data import for EHA ###

#remove(list = ls())

library(igraph); library(magrittr); library(dplyr); library(ggplot2); require(RCurl); library(readr)

AssocsBase <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% data.frame()
HostTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/hosts.csv") %>% data.frame()
VirusTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/viruses.csv") %>% data.frame()

names(AssocsBase)[1:2] <- c("Virus", "Host")
AssocsBase <- mutate(AssocsBase, Virus = as.factor(Virus), Host = as.factor(Host))

AssocsBase2 <- AssocsBase
AssocsBase2 <- droplevels(AssocsBase[!AssocsBase$Host == "Homo_sapiens"&
  !AssocsBase$Virus == "Rabies_virus",])

# Making bipartite projections ####

AssocsTraits <- AssocsBase2[,1:2]

m <- table(AssocsTraits)
M <- as.matrix(m)

bipgraph <- graph.incidence(M, weighted = T)

Virusgraph <- bipartite.projection(bipgraph)$proj1
Hostgraph <- bipartite.projection(bipgraph)$proj2

VirusAdj <- as.matrix(get.adjacency(Virusgraph, attr = "weight"))
diag(VirusAdj) <- table(AssocsBase2$Virus)
VirusA <- matrix(rep(table(AssocsBase2$Virus), nrow(VirusAdj)), nrow(VirusAdj))
VirusB <- matrix(rep(table(AssocsBase2$Virus), each = nrow(VirusAdj)), nrow(VirusAdj))
VirusAdj2 <- VirusAdj/(VirusA + VirusB - VirusAdj)
VirusAdj3 <- VirusAdj/(VirusA) # Asymmetrical

HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(AssocsBase2$Host)
HostA <- matrix(rep(table(AssocsBase2$Host), nrow(HostAdj)), nrow(HostAdj))
HostB <- matrix(rep(table(AssocsBase2$Host), each = nrow(HostAdj)), nrow(HostAdj))
HostAdj2 <- HostAdj/(HostA + HostB - HostAdj)
HostAdj3 <- HostAdj/(HostA)

# Deriving metrics from the networks ####

Hosts <- data.frame(Sp = names(V(Hostgraph)),
                    Degree = degree(Hostgraph),
                    Eigenvector = eigen_centrality(Hostgraph)$vector,
                    Kcore = coreness(Hostgraph),
                    Between = betweenness(Hostgraph))

Viruses <- data.frame(Sp = names(V(Virusgraph)),
                      Degree = degree(Virusgraph),
                      Eigenvector = eigen_centrality(Virusgraph)$vector,
                      Kcore = coreness(Virusgraph),
                      Between = betweenness(Virusgraph))

Hosts <- merge(Hosts, HostTraits, by.x = "Sp", by.y = "hHostNameFinal", all.x = T)
Viruses <- merge(Viruses, VirusTraits, by.x = "Sp", by.y = "vVirusNameCorrected", all.x = T)

Hosts <- Hosts %>% dplyr::rename(hDom = hWildDomFAO)

Domestics <- Hosts[Hosts$hDom == "domestic", "Sp"]
Wildlife <- Hosts[Hosts$hDom == "wild", "Sp"]

DomesticViruses <- as.factor(AssocsBase[AssocsBase$Host %in% Domestics, "Virus"])
WildlifeViruses <- as.factor(AssocsBase[AssocsBase$Host %in% Wildlife, "Virus"])
HumanViruses <- as.factor(AssocsBase[AssocsBase$Host == "Homo_sapiens", "Virus"])

ZoonoticViruses <- intersect(HumanViruses, WildlifeViruses)

AssocsTraits <- merge(AssocsTraits, HostTraits, by.x = "Host", by.y = "hHostNameFinal", all.x = T)
AssocsTraits <- merge(AssocsTraits, VirusTraits, by.x = "Virus", by.y = "vVirusNameCorrected", all.x = T)

AssocsTraits$Domestic <- ifelse(AssocsTraits$Host%in%Domestics,1,0)
AssocsTraits$Wildlife <- ifelse(AssocsTraits$Host%in%Wildlife,1,0)
AssocsTraits$Zoonosis <- ifelse(AssocsTraits$Virus%in%ZoonoticViruses,1,0)

Hosts <- Hosts %>% 
  mutate(
    Domestic = ifelse(Sp %in% Domestics, 1, 0),
    Wildlife = ifelse(Sp %in% Wildlife, 1, 0),
    hZoonosisCount = c(table(AssocsTraits[AssocsTraits$Virus%in%ZoonoticViruses,"Host"])),
    Records = c(table(AssocsTraits$Host))
  ) %>%
  mutate(hZoonosisProp = hZoonosisCount/Records)

Viruses <- Viruses %>%
  mutate(
    Human = case_when(
      Sp %in% HumanViruses ~ 1,
      TRUE ~ 0),
    
    Domestic = case_when(
      Sp %in% DomesticViruses ~ 1,
      TRUE ~ 0),
    
    Wildlife = case_when(
      Sp %in% WildlifeViruses ~ 1,
      TRUE ~ 0),
    
    DomesticCount = c(table(AssocsTraits[AssocsTraits$Domestic == 1,"Virus"])),
    WildlifeCount = c(table(AssocsTraits[AssocsTraits$Wildlife == 1,"Virus"])),
    ZoonosisCount = c(table(AssocsTraits[AssocsTraits$Zoonosis == 1,"Virus"])),
    HumanCount = c(table(AssocsBase[AssocsBase$Host == "Homo_sapiens", "Virus"]))[as.character(Viruses$Sp)],
    
    Records = c(table(AssocsTraits$Virus))
  ) %>% mutate(
    
    
    PropDomestic = DomesticCount/Records,
    PropWildlife = WildlifeCount/Records,
    PropHuman = HumanCount/Records,
    
    PropZoonosis = ZoonosisCount/Records
    
  )


Viruses$HumDomWild <- factor(with(Viruses, 
                                  paste(ifelse(Human,"Human",""), 
                                        ifelse(Domestic,"Domestic",""), 
                                        ifelse(Wildlife,"Wild",""), sep = "")))

vCentrality <- c("vDegree", "vEigenvector", "vCore")
vDists <- c("gaussian", "beta", "binomial")

Viruses$vDegree <- kader:::cuberoot(Viruses$Degree)
Viruses$vEigenvector <- kader:::cuberoot(Viruses$Eigenvector)
Viruses$vCore <- ifelse(Viruses$Kcore==max(Viruses$Kcore),1,0)

#Viruses$vEigenvector[round(Viruses$vEigenvector,4)==0] <- 0.0001 # Needed for a beta-distribution
#Viruses$vEigenvector[round(Viruses$vEigenvector,4)==1] <- 0.999

Viruses$vGenomeAveLengthLn <- log(Viruses$vGenomeAveLength)
Viruses$vPubMedCitesLn <- log(Viruses$vPubMedCites + 1)

# Loading functions, determining themes ####

#devtools::install_github("gfalbery/ggregplot")
library(ggregplot); library(ggplot2); library(RColorBrewer)

ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

AlberTheme <- theme_bw() +
  theme(axis.title.x = element_text(vjust = -0.35, 
                                    size = 12, 
                                    colour = "black"), 
        axis.title.y = element_text(vjust = 1.2, 
                                    size = 12, 
                                    colour = "black"),
        strip.background = element_rect(fill = "white", colour = "dark grey"))

theme_set(AlberTheme)

