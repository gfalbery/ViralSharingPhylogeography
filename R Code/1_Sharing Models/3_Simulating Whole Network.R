
# Simulating on the full network ####

# Rscript "R Code/1_Sharing Models/3_Simulating Whole Network.R"

#source("R Code/00_Master Code.R")

library(MCMCglmm); library(tidyverse); library(Matrix); library(parallel); library(mgcv)

N = nrow(AllMammaldf)

load("~/Albersnet/Output Files/BAMList.Rdata")
load("~/Albersnet/Output Files/BAMList.Rdata")

SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
SpCoef <- BAMList[[1]]$coef[SpCoefNames]

Divisions = round(seq(0, nrow(AllMammaldf), length = 21))

print("Prediction Effects!")

if(file.exists("Output Files/AllPredictions1b.Rdata")) load("Output Files/AllPredictions1b.Rdata") else{
  
  FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")
  
  AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)
  AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), log(FinalHostMatrix$MinCites.Sp2+1)))
  AllMammaldf$Domestic <- 0
  
  AllPredictions1b <- predict.bam(BAMList[[1]], 
                                  newdata = AllMammaldf, # %>% select(-Spp),
                                  type = "terms",
                                  exclude = "Spp")
  
  AllPredictions1b <- mclapply(2:21, function(i){
    predict.bam(BAMList[[1]], 
                newdata = AllMammaldf[(Divisions[i-1]+1):Divisions[i],], 
                type = "terms",
                exclude = "Spp")
  })
  
  save(AllPredictions1b, file = "Output Files/AllPredictions1b.Rdata")
  
}

print("Predicting All Links!")

if(file.exists("Output Files/AllPredList.Rdata")) load("Output Files/AllPredList.Rdata") else{
  
  AllIntercept <- attr(AllPredictions1b[[1]], "constant")
  
  AllPredictions <- lapply(AllPredictions1b, as.data.frame) %>% bind_rows
  
  AllPredictions[,"Intercept"] <- AllIntercept
  
  AllPredList <- parallel::mclapply(1:100, function(x){ # to do something non-specific
    
    AllPredictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
      sample(SpCoef, N, replace = T)
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(AllPredictions)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  save(AllPredList, file = "Output Files/AllPredList.Rdata")
  
}

PredDF1 <- data.frame(AllPredList)

print("Simulating All Networks!")

# Simulating the network #####

if(file.exists("Output Files/AllSims.Rdata")) load("Output Files/AllSims.Rdata") else{
  
  AllSims <- parallel::mclapply(1:length(AllPredList), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = length(AllMammals), 
                     ncol = length(AllMammals))
    
    AssMat[lower.tri(AssMat)] <- round(AllPredList[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                             union(AllMammaldf$Sp,AllMammaldf$Sp2))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  if(length(which(sapply(AllSims, is.null)))>0){
    AllSims <- AllSims[-which(sapply(AllSims, is.null))]
    print("Something went wrong UGH")
    print(paste("New Length = ", length(AllSims)))
  }
  
  save(AllSims, file = "Output Files/AllSims.Rdata")
  
}

# Making summed matrix ####

print("Summing Matrix!")

if(file.exists("Output Files/AllSums.Rdata")) load("Output Files/AllSums.Rdata") else{
  
  AllPredDF <- AllPredList %>% as.data.frame()
  
  #AllPredSums <- apply(AllPredDF,1,sum)
  AllPredSums <- logistic(rowSums(AllPredictions))
  
  AssMat <- matrix(NA, 
                   nrow = length(AllMammals), #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(AllMammals)) #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[lower.tri(AssMat)] <- AllPredSums
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  
  dimnames(AssMat) <- list(AllMammals,
                           AllMammals)
  
  AllSums <- as(AssMat, "dgCMatrix")
  
  save(AllSums, file = "Output Files/AllSums.Rdata")
  
}

# Making into Graphs ####

print("Making All Graphs!")

if(file.exists("Output Files/AllSimGs.Rdata")) load("Output Files/AllSimGs.Rdata") else{
  
  AllSimGs <- parallel::mclapply(1:length(AllSims), function(i){
    
    graph.adjacency(AllSims[[i]], mode = "undirected", diag = F)
    
  }, mc.cores = 10)
  
  if(length(which(sapply(AllSims, is.null)))>0){
    AllSimGs <- AllSimGs[-which(sapply(AllSims, is.null))]
    print("Something went wrong UGH")
  }
  
  save(AllSimGs, file = "Output Files/AllSimGs.Rdata")
  
}

print(Sys.time())