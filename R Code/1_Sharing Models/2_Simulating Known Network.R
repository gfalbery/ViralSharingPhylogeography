
# STAN Model Output ####

# Rscript "R Code/1_Sharing Models/2_Simulating Known Network.R"

library(ggregplot); library(parallel); library(igraph);
library(mgcv); library(tidyverse)

if(file.exists("Output Files/KnownSimGraphs.Rdata")) load("Output Files/KnownSimGraphs.Rdata") else{
  
  if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")
  
  load("Output Files/BAMList.Rdata")
  
  Resps <- c("VirusBinary")#,"RNA","DNA","Vector","NVector")
  
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
  FinalHostMatrix$PredVirus1Q <- cut(FinalHostMatrix$PredVirus1,
                                     breaks = c(-1:10/10),
                                     labels = c(0:10/10))
  
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
  
  Predictions <- Predictions1b %>% as.data.frame
  
  N = nrow(FinalHostMatrix)
  
  PredList1b <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    Predictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
      sample(SpCoef, N, replace = T)
    
    Predictions[,"Intercept"] <- Intercept1b
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1b <- data.frame(PredList1b)
  FinalHostMatrix$PredVirus1b <- apply(PredDF1b, 1, mean)
  FinalHostMatrix$PredVirus1bQ <- cut(FinalHostMatrix$PredVirus1b,
                                      breaks = c(-1:10/10),
                                      labels = c(0:10/10))
  
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
  
  Predictions <- Predictions1c %>% as.data.frame
  
  N = nrow(FinalHostMatrix)
  
  PredList1c <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    Predictions[,"Intercept"] <- Intercept1c
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1c <- data.frame(PredList1c)
  FinalHostMatrix$PredVirus1c <- logistic(rowSums(Predictions)) # apply(PredDF1b, 1, mean)
  FinalHostMatrix$PredVirus1cQ <- cut(FinalHostMatrix$PredVirus1c,
                                      breaks = c(-1:10/10),
                                      labels = c(0:10/10))
  
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
