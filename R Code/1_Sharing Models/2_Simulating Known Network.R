
# STAN Model Output ####

# Rscript "R Code/1_Sharing Models/2_Simulating Known Network.R"

library(rstan); library(reskew); library(ggregplot); library(parallel); library(igraph);
library(mgcv); library(tidyverse)

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

load("Output Files/BAMList.Rdata")

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

Resps <- c("VirusBinary")#,"RNA","DNA","Vector","NVector")

DataList <- list()

for(r in 1:length(Resps)){
  
  print(Resps[r])
  
  DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r]))
  
  DataList[[Resps[r]]]$Sp <- factor(DataList[[Resps[r]]]$Sp, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  DataList[[Resps[r]]]$Sp2 <- factor(DataList[[Resps[r]]]$Sp2, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  
  MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  
  SppMatrix = MZ1 + MZ2
  
  DataList[[Resps[[r]]]]$Spp <- SppMatrix
  DataList[[Resps[[r]]]]$Cites <- rowSums(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
  DataList[[Resps[[r]]]]$MinCites <- apply(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
  DataList[[Resps[[r]]]]$Domestic <- ifelse(rowSums(cbind(2- FinalHostMatrix$hDom %>% as.factor %>% as.numeric,
                                                          2- FinalHostMatrix$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)
}

# Doing the simulating with random effects #####

PredList1 <- list()

Predictions1 <- predict.gam(BAMList[[1]], 
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
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[lower.tri(AssMat)] <- round(PredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

SimGraphs1 <- mclapply(1:length(PredList1), function(i){
  
  graph.adjacency(SimNets1[[i]], mode = "undirected")
  
}, mc.cores = 10)

# Doing the simulating with random effects #####

PredList1b <- list()

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

save(SimGraphs1, SimGraphs1b, file = "Output Files/KnownSimGraphs.Rdata")

