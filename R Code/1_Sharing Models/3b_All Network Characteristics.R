
# 3b_Getting Whole Network Characteristics ####

# Rscript "R Code/1_Sharing Models/3b_All Network Characteristics.R"

library(igraph); library(tidyverse); library(ggregplot); library(parallel); library(SpRanger); library(Matrix)

#source("R Code/00_Master Code.R")

print("Getting All Network Characteristics")

load("Output Files/AllSimGs.Rdata")
load("Output Files/AllSims.Rdata")

print("Deriving Characteristics")

if(file.exists("Output Files/AllNetworkStats.Rdata")) load("Output Files/AllNetworkStats.Rdata") else{
  
  AllPredNetwork <- mclapply(AllSimGs, AllNetworkStats, mc.cores = 10)
  
  save(AllPredNetwork, file = "Output Files/AllNetworkStats.Rdata")
  
}

print("Getting links")

AllDegrees <- map(AllPredNetwork, "Degree") %>% bind_cols()

AllPredDegree = apply(AllDegrees, 1, mean)

AllPredDegrees <- data.frame(AllPredDegree, Sp = AllMammals)

Hosts <- left_join(Hosts, AllPredDegrees, by = "Sp")
# GGally::ggpairs(Hosts %>% select(contains("Degree")), lower = list(continuous = "smooth"))

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 <- left_join(Panth1, AllPredDegrees, by = "Sp")

load("data/FullPolygons.Rdata")

FullPolygons <- left_join(FullPolygons, AllPredDegrees, by = c("Host" = "Sp"), all.x = T)

print("Dividing within- and between-order links")

# Identifying within-order links ####

Panth1 <- Panth1 %>% slice(order(Panth1$Sp)) %>% 
  filter(Sp%in%V(AllSimGs[[1]])$name) %>% droplevels

hOrderList <- mclapply(1:length(AllSimGs), function(i){
  
  a <- AllSimGs[[i]]
  
  lapply(levels(Panth1$hOrder), function(x){
    
    distances(a, v = V(a)[Panth1$hOrder == as.character(x)], 
              to = V(a)[Panth1$hOrder == as.character(x)])
    
  })
  
}, mc.cores = 10) #%>% bind_rows()

InDegrees <- lapply(hOrderList, function(a) lapply(a, function(b) rowSums(b==1)))

InDegrees2 <- lapply(InDegrees, unlist) %>% bind_cols %>% apply(1, mean)

InDegrees2 <- data.frame(Sp = lapply(InDegrees, unlist)[[1]] %>% names,
                         InDegree = InDegrees2)

Panth1 <- Panth1 %>% left_join(InDegrees2)

Panth1$OutDegree <- with(Panth1, AllPredDegree-InDegree)

EIDSpecies <- read.csv("data/EID/SpeciesInteractions_EID2.csv", header = T)

#EIDSpecies #<- EIDSpecies %>% filter(Sequences.count>0)

EIDSpecies <- EIDSpecies %>% mutate(Carrier = str_replace(Carrier, " ", "_"))

substr(EIDSpecies$Carrier,1,1) = toupper(substr(EIDSpecies$Carrier,1,1))

Panth1 <- Panth1 %>% mutate(Obs = ifelse(Sp %in% FHN, 1, 0),
                            ZoonoticHost = ifelse(Sp %in% Hosts[Hosts$hZoonosisCount>0,"Sp"], 1, 0)) %>%
  mutate(EIDObs = ifelse(Sp %in% EIDSpecies$Carrier, 1, 0)) %>%
  mutate(JustEID = ifelse(EIDObs==1&Obs==0, 1, 0)) 

save(Panth1, file = "Output Files/Panth1.Rdata")


scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

