
# Simulating on the full network ####

library(MCMCglmm); library(tidyverse); library(Matrix)

invlogit <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeOverlap))
AllMammals <- AllMammals[order(AllMammals)]

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj1[AllMammals,AllMammals]),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
)

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

MX1 = model.matrix( ~ Space + Phylo + Space:Phylo, data = AllMammaldf) %>%
  as.matrix %>% as("dgCMatrix")

# Trying it without random effects ####

i = 2

AllPredList <- list()

ClusterMCMC <- BinModelList[1:10 + 10*(i-1)] %>% 
  lapply(function(a) as.data.frame(a$Sol)) %>% 
  bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

for(x in 1:1000){
  if(x%%10==0) print(x)
  RowSampled <- RowsSampled[x]
  FXSample <- ClusterMCMC[RowSampled, c("(Intercept)","Space","Phylo2","Space:Phylo2")]
  Output <- c(FXSample %*% t(MX1))
  PZero <- rbinom(length(Output[[1]]@x), 1, invlogit(Output[[1]]@x))
  AllPredList[[x]] <- PZero
}

AllPredDF <- as.data.frame(AllPredList)
AllMammaldf$PredVirus <- apply(AllPredDF,1, function(a) a %>% mean)

FinalHostMatrix$PredVirusQ <- cut(AllMammaldf$PredVirus,
                                  breaks = c(-1:10/10),
                                  labels = c(0:10/10))

# Simulating the network #####

AllSimNets <- AllSimGraphs <- list()

for(i in 1:length(AllPredList)){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList2[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  AllSimNets[[i]] <- as(AssMat, "dgCMatrix")
  
  AllSimGraphs[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
  if(i%%10==0) print(i)
  
}

AllDegdf <- sapply(AllSimGraphs, function(a) degree(a)) %>% as.data.frame
AllEigendf <- sapply(AllSimGraphs, function(a) eigen_centrality(a)$vector) %>% as.data.frame

AllPredDegrees <- apply(AllDegdf, 1, mean)
AllPredEigen <- apply(AllEigendf, 1, mean)

