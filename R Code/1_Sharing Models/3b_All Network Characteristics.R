
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

Sp = rownames(AllSims[[1]])

AllPredDegrees <- data.frame(AllPredDegree, Sp)

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

for(z in c("Method1", "Method2")){
  
  if(z == "Method1"){
    
    hOrderList <- mclapply(1:length(AllSimGs), function(i){
      lapply(levels(Panth1$hOrder), function(x) induced_subgraph(AllSimGs[[i]], 
                                                                 as.character(Panth1$hOrder) == x))
    }, mc.cores = 10)
    
    InDegree <- lapply(hOrderList, function(a) lapply(a, function(b) degree(b)) %>% unlist) %>% bind_cols
    InDegrees <- apply(InDegree, 1, mean)
    InNames <- lapply(hOrderList[[1]], function(a) V(a)$name) %>% unlist
    InDegDF <- data.frame(InDegree = InDegrees,
                          Sp = InNames)
    
    Panth1 <- left_join(Panth1, InDegDF, by = "Sp")
    
  } else {
    
    hOrderList <- mclapply(levels(Panth1$hOrder), function(i){
      
      map(AllSimGs, function(j){
        
        distances(j, v = V(j)[Panth1$hOrder == as.character(i)], 
                  to = V(j)[Panth1$hOrder == as.character(i)])
        
      }) %>% bind_cols()
      
    }, mc.cores = 10)
    
    hOrderList <- mclapply(1:length(AllSimGs), function(i){
      
      a <- AllSimGs[[i]]
      
      lapply(levels(Panth1$hOrder), function(x){
        
        distances(a, v = V(a)[Panth1$hOrder == as.character(x)], 
                  to = V(a)[Panth1$hOrder == as.character(x)])
        
      })
      
    }, mc.cores = 10) %>% bind_rows()
    
    InDegrees <- lapply(hOrderList, function(a) lapply(a, function(b) rowSums(b==1)))
    
    InDegrees2 <- lapply(1:nlevels(Panth1$hOrder), 
                         function(a) map(InDegrees, a)) #%>% bind_rows()
    
  }
  
  
}

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


Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>% right_join(Panth1, by = "Sp") %>%
  mutate(CentredDegree = )




