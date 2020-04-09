
# Rscript "4_ Whole Network Characteristics.R" ####

# Rscript "R Code/1_Sharing Models/3b_All Network Characteristics.R"

library(igraph); library(tidyverse); library(ggregplot); library(parallel); library(SpRanger); library(Matrix)

print("Getting All Network Characteristics")

load("~/Albersnet/Output Files/AllSimGs.Rdata")
load("~/Albersnet/Output Files/AllSums.Rdata")

if(file.exists("Output Files/AllNetworkStats.Rdata")) load("Output Files/AllNetworkStats.Rdata") else{
  
  AllPredNetwork <- mclapply(AllSimGs, AllNetworkStats, mc.cores = 10)
  
  save(AllPredNetwork, file = "Output Files/AllNetworkStats.Rdata")
  
}

print("Getting links")

AllDegrees <- map(AllPredNetwork, "Degree") %>% bind_cols()
AllPredDegree <- apply(AllDegrees, 1, mean)
AllPredDegrees <- data.frame(AllPredDegree, Sp = AllMammals)

Hosts <- left_join(Hosts, AllPredDegrees, by = "Sp")

save(Hosts, file = "Output Files/Hosts.Rdata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 <- left_join(Panth1, AllPredDegrees, by = "Sp")

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
EIDSpecies <- EIDSpecies %>% mutate(Carrier = str_replace(Carrier, " ", "_") %>% CamelConvert,
                                    Cargo = str_replace(Cargo, " ", "_"))

substr(EIDSpecies$Carrier,1,1) <- toupper(substr(EIDSpecies$Carrier,1,1))
substr(EIDSpecies$Cargo,1,1) <- toupper(substr(EIDSpecies$Cargo,1,1))

Panth1 <- Panth1 %>% mutate(Obs = ifelse(Sp %in% FHN, 1, 0),
                            ZoonoticHost = ifelse(Sp %in% Hosts[Hosts$hZoonosisCount>0,"Sp"], 1, 0)) %>%
  mutate(EIDObs = ifelse(Sp %in% EIDSpecies$Carrier, 1, 0)) %>%
  mutate(JustEID = ifelse(EIDObs==1&Obs==0, 1, 0)) %>%
  mutate(Subset = paste0(Obs, EIDObs))

save(Panth1, file = "Output Files/Panth1.Rdata")

# Making EID sharing network ####

AssocsEID <- EIDSpecies[EIDSpecies$Cargo.classification == "Virus"&
                          EIDSpecies$Carrier%in%AllMammals,
                        c("Cargo", "Carrier")] %>%
  droplevels %>%
  slice(order(Carrier))

m <- table(AssocsEID)
M <- as.matrix(m)

bipgraph <- graph.incidence(M, weighted = NULL)

EIDgraph <- bipartite.projection(bipgraph)$proj2

EIDAdj <- as.matrix(get.adjacency(EIDgraph, attr = "weight"))
diag(EIDAdj) <- table(AssocsEID$Carrier)
Remove <- which(rowSums(EIDAdj)==diag(EIDAdj))
EIDAdj <- EIDAdj[-Remove,-Remove]

EIDAdj[EIDAdj>0] <- 1

EIDMammals <- intersect(AllMammals, rownames(EIDAdj))

EIDCordf <- data.frame(
  EIDConnected = c(EIDAdj[EIDMammals,EIDMammals][lower.tri(EIDAdj[EIDMammals,EIDMammals])]),
  PredNetwork = c(AllSums[EIDMammals,EIDMammals][lower.tri(EIDAdj[EIDMammals,EIDMammals])])
)

EIDCordf[,c("Sp","Sp2")] <- expand.grid(EIDMammals, EIDMammals)[lower.tri(EIDAdj[EIDMammals,EIDMammals])]

# Summarising predicted degree by grid square ####

library(tidyverse); library(raster); library(colorspace)

if(file.exists("Output Files/GridDegree.Rdata")) load("Output Files/GridDegree.Rdata") else{
  
  load("~/LargeFiles/MammalStackFullMercator.Rdata")
  
  RasterListb <- lapply(1:length(AllMammals), function(a){
    
    if(a %% 500==0) print(a)
    
    MammalStackFull[[AllMammals[a]]]
    
  })
  
  ToSkip <- which(sapply(RasterListb, is.null))
  
  names(RasterListb) <- AllMammals
  
  OverlapSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  
  for(i in 1:length(RasterListb)){
    
    if(i %% 500 == 0) print(i)
    SubSums <- raster::getValues(RasterListb[[i]]) > 0 %>% as.numeric
    SubSums[is.na(SubSums)] <- 0
    OverlapSums <- OverlapSums + SubSums
    
  }
  
  AllDegree <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  InDegree <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  OutDegree <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  
  for(i in 1:length(RasterListb)){
    
    if(i %% 500 == 0) print(i)
    SubSums <- raster::getValues(RasterListb[[i]]) > 0 %>% as.numeric
    SubSums[is.na(SubSums)] <- 0
    AllDegree <- AllDegree + SubSums*Panth1[Panth1$Sp==AllMammals[i],"AllPredDegree"]
    InDegree <- InDegree + SubSums*Panth1[Panth1$Sp==AllMammals[i],"InDegree"]
    OutDegree <- OutDegree + SubSums*Panth1[Panth1$Sp==AllMammals[i],"OutDegree"]
    
  }
  
  DegreeDF <- data.frame(
    X = rep(1:ncol(MammalStackFull[[1]]), nrow(MammalStackFull[[1]])),
    Y = rep(nrow(MammalStackFull[[1]]):1, each = ncol(MammalStackFull[[1]])),
    
    Richness = OverlapSums,
    AllDegree = AllDegree,
    InDegree = InDegree,
    OutDegree = OutDegree
    
  ) 
  
  UniversalBlank <- raster::raster("Iceberg Input Files/UniversalBlank.tif")
  Land = which(raster::values(UniversalBlank)==0)
  Sea = which(is.na(raster::values(UniversalBlank)))
  
  DegreeDF <- DegreeDF[-Sea,]
  
  GridDegree <- DegreeDF %>% 
    mutate_at(vars(contains("Degree")), function(a) a/DegreeDF$Richness) %>% 
    mutate_at(vars(contains("Degree")), function(a) ifelse(a==0|is.na(a), min(a[a>0], na.rm = T), a)) %>% 
    tidyr::gather(key = "Metric", value = "Degree", contains("Degree"))
  
  save(GridDegree, file = "Output Files/GridDegree.Rdata")
  
}

# Dissecting between-order links ####

Combos <- t(combn(levels(Panth1$hOrder),2)) %>% as.data.frame() %>%
  rename(Order1 = V1, Order2 = V2)

hComboList <- list()

for(i in 1:nrow(Combos)){
  
  print(i)
  
  a = Combos[i,] %>% unlist
  
  b1 = Panth1[Panth1$hOrder%in%a[1],"Sp"] %>% intersect(AllMammals)
  b2 = Panth1[Panth1$hOrder%in%a[2],"Sp"] %>% intersect(AllMammals)
  
  if(length(b1)>0&length(b2)>0){
    
    FocalNet <- AllSums[b1,b2] %>% as.matrix
    
    hComboList[[paste(Combos[i,1],Combos[i,2], sep = ".")]] <- data.frame(Degree = c(rowSums(FocalNet), colSums(FocalNet)),
                                  Iteration = i,
                                  Sp = c(b1, b2),
                                  Order = c(rep(a[1], length(b1)), rep(a[2], length(b2))),
                                  Group = rep(1:2, c(length(b1),length(b2))))
  }
}

hComboList <- hComboList %>% bind_rows(.id = "Combo")

hUniteList <- list()

for(i in levels(Panth1$hOrder)){
  
  print(i)
  
  b1 = Panth1[Panth1$hOrder == i, "Sp"] %>% intersect(AllMammals)
  
  if(length(b1)>0){
    FocalNet <- AllSums[b1,b1] %>% as.matrix
    
    hUniteList[[i]] <- data.frame(Degree = c(rowSums(FocalNet)),
                                  Iteration = i,
                                  Sp = b1,
                                  Order = i)
  }
}

hUniteList <- hUniteList %>% bind_rows()

OutDegrees <- hComboList %>% group_by(Sp) %>% summarise(OutDegree = sum(Degree))
InDegrees <- hUniteList %>% group_by(Sp) %>% summarise(InDegree = sum(Degree))
AllDegrees <- OutDegrees %>% left_join(InDegrees) %>%
  mutate(AllPredDegree = OutDegree + InDegree)

Panth1 <- Panth1 %>% #dplyr::select(-c("AllPredDegree", "InDegree", "OutDegree")) %>%
  left_join(AllDegrees, by = "Sp")

OrderLevelLinks <- Panth1 %>% #dplyr::select(-c("AllPredDegree", "InDegree", "OutDegree")) %>%
  left_join(AllDegrees, by = "Sp") %>%
  group_by(hOrder) %>%
  summarise(HostNumber = n(),
            OutDegree = mean(OutDegree, na.rm = T),
            InDegree = mean(InDegree, na.rm = T),
            AllPredDegree = mean(AllPredDegree, na.rm = T)) %>%
  tidyr::gather(key = "Metric", value = "Degree", contains("Degree"))


