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

bipgraph <- graph.incidence(M, weighted = NULL)

Hostgraph <- bipartite.projection(bipgraph)$proj2

HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(AssocsBase2$Host)
Remove <- which(rowSums(HostAdj)==diag(HostAdj))
HostAdj <- HostAdj[-Remove,-Remove]

# Deriving metrics from the networks ####

Hosts <- data.frame(Sp = names(V(Hostgraph)),
                    Degree = degree(Hostgraph),
                    Eigenvector = eigen_centrality(Hostgraph)$vector,
                    Kcore = coreness(Hostgraph),
                    Between = betweenness(Hostgraph))

Hosts <- merge(Hosts, HostTraits, by.x = "Sp", by.y = "hHostNameFinal", all.x = T)
Hosts <- Hosts %>% dplyr::rename(hDom = hWildDomFAO)

Domestics <- Hosts[Hosts$hDom == "domestic", "Sp"]
Wildlife <- Hosts[Hosts$hDom == "wild", "Sp"]

AssocsTraits <- merge(AssocsTraits, HostTraits, by.x = "Host", by.y = "hHostNameFinal", all.x = T)

AssocsTraits$Domestic <- ifelse(AssocsTraits$Host%in%Domestics,1,0)
AssocsTraits$Wildlife <- ifelse(AssocsTraits$Host%in%Wildlife,1,0)

ZoonoticViruses <- AssocsBase %>% filter(Host == "Homo_sapiens") %>% select(Virus)

Hosts <- Hosts %>% 
  mutate(
    Domestic = ifelse(Sp %in% Domestics, 1, 0),
    Wildlife = ifelse(Sp %in% Wildlife, 1, 0),
    hZoonosisCount = c(table(AssocsTraits[AssocsTraits$Virus%in%ZoonoticViruses$Virus,"Host"])),
    Records = c(table(AssocsTraits$Host))
  )

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

