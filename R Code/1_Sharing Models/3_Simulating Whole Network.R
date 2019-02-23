
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

Predictions <- predict(BAMList[[1]], newdata = AllMammaldf)

N = nrow(AllMammaldf)

AllPredList <- parallel::mclapply(1:100, function(x){ # to do something non-specific
  
  BinPred <- rbinom(n = N,
                    prob = logistic(Predictions),
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

AllPrev = sapply(AllPredList, Prev)

AllDegdf <- sapply(AllSimGs, function(a) degree(a))
AllPredDegrees <- apply(AllDegdf, 1, mean)

AllComponents <- sapply(AllSimGs, function(b) components(b)$no)

Hosts[,"AllPredDegree"] <- AllPredDegrees[as.character(Hosts$Sp)]

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")
Panth1[,"AllPredDegree"] <- AllPredDegrees[as.character(Panth1$Sp)]

FullPolygons[,"AllPredDegree"] <- AllPredDegrees[as.character(FullPolygons$Host)]

FullCentroids <-  FullPolygons %>% group_by(Host) %>% 
  summarise(LongMean = (max(long) + min(long))/2,
            LatMean = (max(lat) + min(lat))/2) %>% data.frame

rownames(FullCentroids) <- FullCentroids$Host

Panth1[,c("LongMean","LatMean")] <- FullCentroids[Panth1$Sp,c("LongMean","LatMean")]

# Identifying within-order links ####

Panth1 <- Panth1 %>% slice(order(Panth1$Sp)) %>% 
  filter(Sp%in%V(AllSimGs[[1]])$name) %>% droplevels

hOrderList <- mclapply(1:length(AllSimGs), function(i){
  lapply(levels(Panth1$hOrder), function(x) induced_subgraph(AllSimGs[[i]], 
                                                             as.character(Panth1$hOrder) == x))
})

InDegree = lapply(hOrderList, function(a) lapply(a, function(b) degree(b)) %>% unlist)
InNames <- lapply(hOrderList[[1]], function(a) V(a)$name ) %>% unlist
InDegDF <- data.frame(InDegree)

InDegrees <- apply(InDegDF,1, function(a) mean(a, na.rm = T))

Panth1$InDegree <- InDegrees[Panth1$Sp]
Panth1$OutDegree <- with(Panth1, AllPredDegree-InDegree)

