
# Simulating on the full network ####

library(MCMCglmm); library(tidyverse); library(Matrix); library(parallel)

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj1[AllMammals,AllMammals]),
  Phylo2 = scale(c(tFullSTMatrix[AllMammals,AllMammals])),
  DietSim = c(tFullSTMatrix[AllMammals,AllMammals])
)

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

AllPredList <- list()

N = nrow(AllMammaldf)

load("~/Albersnet/Output Files/BAMList.Rdata")

SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
SpCoef <- BAMList[[1]]$coef[SpCoefNames]

Predictions1b <- predict.bam(BAMList[[1]], 
                             newdata = DataList[[1]], # %>% select(-Spp),
                             type = "terms",
                             exclude = "Spp")

Intercept1b <- attr(Predictions1b, "constant")

Predictions <- Predictions1b %>% as.data.frame

AllPredList <- parallel::mclapply(1:100, function(x){ # to do something non-specific
  
  Predictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
    sample(SpCoef, N, replace = T)
  
  Predictions[,"Intercept"] <- Intercept1b
  
  BinPred <- rbinom(n = N,
                    prob = logistic(rowSums(Predictions)),
                    size  = 1)
  
  BinPred
  
}, mc.cores = 10)

PredDF1 <- data.frame(AllPredList)
AllMammaldf$PredVirus1 <- apply(PredDF1, 1, mean)
AllMammaldf$PredVirus1Q <- cut(AllMammaldf$PredVirus1,
                               breaks = c(-1:10/10),
                               labels = c(0:10/10))

# Simulating the network #####

AllSims <- parallel::mclapply(1:length(AllPredList), function(i){
  
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

AllSimGs <- parallel::mclapply(1:length(AllSims), function(i){
  
  graph.adjacency(AllSims[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

if(length(which(sapply(AllSims, is.null)))>0) AllSimGs <- AllSimGs[-which(sapply(AllSims, is.null))]

save(AllSims, file = "Output Files/AllSims.Rdata")
save(AllSimGs, file = "Output Files/AllSimGs.Rdata")
