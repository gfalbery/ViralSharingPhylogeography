
# STAN Model Output ####

# Rscript "R Code/1_Sharing Models/2_Simulating Known Network.R"

library(ggregplot); library(parallel); library(igraph);
library(mgcv); library(tidyverse)

if(file.exists("Output Files/KnownSimGraphs.Rdata")) load("Output Files/KnownSimGraphs.Rdata") else{
  
  if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")
  
  load("Output Files/BAMList.Rdata")
  
  Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")
  
  # Doing the simulating with random effects #####
  
  PredList1 <- list()
  
  Predictions1 <- predict.bam(BAMList[[1]], 
                              newdata = DataList[[1]],
                              type = "terms")
  
  Intercept1 <- attr(Predictions1, "constant")
  
  Predictions1 <- Predictions1 %>% as.data.frame()
  
  Predictions1$Intercept <- Intercept1
  
  N = nrow(DataList[[1]])
  
  PredList1 <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions1)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1 <- data.frame(PredList1)
  FinalHostMatrix$PredVirus1 <- apply(PredDF1, 1, mean)
  
  SimNets1 <- mclapply(1:length(PredList1), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = nlevels(DataList[[1]]$Sp), 
                     ncol = nlevels(DataList[[1]]$Sp))
    
    AssMat[lower.tri(AssMat)] <- round(PredList1[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(levels(DataList[[1]]$Sp),
                             levels(DataList[[1]]$Sp))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  SimGraphs1 <- mclapply(1:length(PredList1), function(i){
    
    graph.adjacency(SimNets1[[i]], mode = "undirected")
    
  }, mc.cores = 10)
  
  # Doing the simulating without random effects #####
  
  SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[1]]$coef[SpCoefNames]
  
  Predictions1b <- predict.bam(BAMList[[1]], 
                               newdata = DataList[[1]],# %>% select(-Spp),
                               type = "terms",
                               exclude = "Spp")
  
  Intercept1b <- attr(Predictions1b, "constant")
  
  Predictions1b <- Predictions1b %>% as.data.frame
  Predictions1b[,"Intercept"] <- Intercept1b
  
  N = nrow(FinalHostMatrix)
  
  PredList1b <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    Predictions1b[,"Spp"] <- sample(SpCoef, N, replace = T) + 
      sample(SpCoef, N, replace = T)
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions1b)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1b <- data.frame(PredList1b)
  FinalHostMatrix$PredVirus1b <- apply(PredDF1b, 1, mean)
  
  SimNets1b <- mclapply(1:length(PredList1b), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                     ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
    
    AssMat[lower.tri(AssMat)] <- round(PredList1b[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                             union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  SimGraphs1b <- mclapply(1:length(PredList1b), function(i){
    
    graph.adjacency(SimNets1b[[i]], mode = "undirected")
    
  }, mc.cores = 10)
  
  # Doing the simulating with only random effects #####
  
  SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[1]]$coef[SpCoefNames]
  
  Predictions1c <- predict.bam(BAMList[[1]], 
                               newdata = DataList[[1]] %>% 
                                 mutate_at(vars(Space, Phylo, MinCites), function(a) mean(a)) %>%
                                 mutate(Gz = 0),
                               type = "terms")
  
  Intercept1c <- attr(Predictions1c, "constant")
  
  Predictions1c <- Predictions1c %>% as.data.frame
  Predictions1c[,"Intercept"] <- Intercept1c
  
  N = nrow(FinalHostMatrix)
  
  PredList1c <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions1c)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1c <- data.frame(PredList1c)
  FinalHostMatrix$PredVirus1c <- apply(PredDF1b, 1, mean)
  
  SimNets1c <- mclapply(1:length(PredList1c), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                     ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
    
    AssMat[lower.tri(AssMat)] <- round(PredList1c[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                             union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  SimGraphs1c <- mclapply(1:length(PredList1c), function(i){
    
    graph.adjacency(SimNets1c[[i]], mode = "undirected")
    
  }, mc.cores = 10)
  
  save(SimGraphs1, SimGraphs1b, SimGraphs1c, file = "Output Files/KnownSimGraphs.Rdata")
  
  print(Sys.time())
  
}

# Getting network characteristics ####

library(igraph); library(tidyverse); library(ggregplot); library(parallel); library(SpRanger)

if(file.exists("Output Files/KnownNetworkStats.Rdata")) load("Output Files/KnownNetworkStats.Rdata") else{
  
  print("Observed")
  ObsNetwork <- AllNetworkStats(Hostgraph)
  
  print("1!")
  PredNetwork1 <- mclapply(SimGraphs1, AllNetworkStats, mc.cores = 10)
  
  print("1b!")
  PredNetwork1b <- mclapply(SimGraphs1b, AllNetworkStats, mc.cores = 10)
  
  print("1c!")
  PredNetwork1c <- mclapply(SimGraphs1c, AllNetworkStats, mc.cores = 10)
  
  save(ObsNetwork, PredNetwork1, PredNetwork1b, PredNetwork1c, file = "Output Files/KnownNetworkStats.Rdata")
  
}

# Comparison of degree stuff ####

PredDegrees1 <- map(PredNetwork1, "Degree") %>% bind_cols()
PredDegrees1b <- map(PredNetwork1b, "Degree") %>% bind_cols()
PredDegrees1c <- map(PredNetwork1c, "Degree") %>% bind_cols()

PredDegrees <- data.frame(PredDegree1 = apply(PredDegrees1, 1, mean),
                          PredDegree1b = apply(PredDegrees1b, 1, mean),
                          PredDegree1c = apply(PredDegrees1c, 1, mean),
                          Sp = names(PredNetwork1c[[1]]$Degree))

Hosts <- Hosts %>%
  left_join(PredDegrees, by = "Sp")
