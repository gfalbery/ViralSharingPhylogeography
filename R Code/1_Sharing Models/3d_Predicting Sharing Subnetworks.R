
# Generating subnetworks

# Rscript "R Code/1_Sharing Models/3d_Predicting Sharing Subnetworks.R"

source("R Code/00_Master Code.R")

library(MCMCglmm); library(tidyverse); library(Matrix); library(parallel); library(mgcv)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)
tVD <- 1 - VD

AllMammals <- reduce(list(colnames(FullSTMatrix),
                          colnames(FullRangeAdj1),
                          colnames(VD)), 
                     intersect)

AllMammals <- AllMammals[order(AllMammals)]

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj1[AllMammals,AllMammals]),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals]),
  DietSim = c(tVD[AllMammals,AllMammals])
)

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

N = nrow(AllMammaldf)

Divisions = round(seq(0, nrow(AllMammaldf), length = 21))

AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), log(FinalHostMatrix$MinCites.Sp2+1)))
AllMammaldf$Domestic <- 0

load("~/Albersnet/Output Files/BAMList.Rdata")

print("Prediction Effects!")

SubSums <- lapply(2:5, function(r){
  
  print(Resps[r])
  
  SpCoefNames <- names(BAMList[[r]]$coef)[substr(names(BAMList[[r]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[r]]$coef[SpCoefNames]
  
  FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")
  
  AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)
  
  AllPredictions1b <- mclapply(2:21, function(i){
    predict.bam(BAMList[[r]], 
                newdata = AllMammaldf[(Divisions[i-1]+1):Divisions[i],], 
                type = "terms",
                exclude = "Spp")
  }, mc.cores = 10)
  
  print("Predicting All Links!")
  
  AllIntercept <- attr(AllPredictions1b[[1]], "constant")
  
  AllPredictions <- lapply(AllPredictions1b, as.data.frame) %>% bind_rows
  
  AllPredictions[,"Intercept"] <- AllIntercept
  
  AssMat <- matrix(NA, 
                   nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[lower.tri(AssMat)] <- round(logistic(rowSums(AllPredictions)))
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                           union(AllMammaldf$Sp,AllMammaldf$Sp2))
  
  
  return(AssMat)
  
})

names(SubSums) <- Resps[2:5]

save(SubSums, file = "Output Files/SubSums.Rdata")

stop()

AllPredList <- parallel::mclapply(1:100, function(x){ # to do something non-specific
  
  #AllPredictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
  #  sample(SpCoef, N, replace = T)
  
  BinPred <- rbinom(n = N,
                    prob = logistic(rowSums(AllPredictions)),
                    size  = 1)
  
  BinPred
  
}, mc.cores = 10)



PredDF1 <- data.frame(AllPredList)

print("Simulating All Networks!")

# Simulating the network #####


SubSims <- parallel::mclapply(1:length(AllPredList), function(i){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[lower.tri(AssMat)] <- round(AllPredList[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                           union(AllMammaldf$Sp,AllMammaldf$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

if(length(which(sapply(SubSims, is.null)))>0){
  SubSims <- SubSims[-which(sapply(SubSims, is.null))]
  print("Something went wrong UGH")
  print(paste("New Length = ", length(SubSims)))
}



# Making summed matrix ####

print("Summing Matrix!")

AllPredDF <- AllPredList %>% as.data.frame()

AllPredSums <- apply(AllPredDF,1,sum)

AssMat <- matrix(NA, 
                 nrow = 4276, #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                 ncol = 4276) #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))

AssMat[lower.tri(AssMat)] <- AllPredSums
AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
diag(AssMat) <- 0

dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                         union(AllMammaldf$Sp,AllMammaldf$Sp2))

AllSums <- as(AssMat, "dgCMatrix")


save(SubSums, file = "Output Files/SubSimGs.Rdata")


# Making into Graphs ####

print("Making All Graphs!")

SubSimGs <- parallel::mclapply(1:length(SubSims), function(i){
  
  graph.adjacency(SubSims[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

if(length(which(sapply(SubSims, is.null)))>0){
  SubSimGs <- SubSimGs[-which(sapply(SubSims, is.null))]
  print("Something went wrong UGH")
}

save(SubSimGs, file = "Output Files/SubSimGs.Rdata")


print(Sys.time())