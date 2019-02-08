
# Simulating on the full network ####

library(MCMCglmm); library(tidyverse); library(Matrix)

invlogit <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj1[AllMammals,AllMammals]),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
)

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

# Trying it without random effects ####

AllPredList <- list()


MX1 = model.matrix( ~ Space + Phylo + Space:Phylo, data = AllMammaldf) %>%
  as.matrix %>% as("dgCMatrix")

N = nrow(AllMammaldf)

ClusterMCMC <- p$df

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

for(x in 1:length(RowsSampled)){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  XFX <- p$df[RowSampled, XBetas] %>% unlist
  XPredictions <- c(XFX %*% t(MX1))
  
  ZPredictions2a <- rnorm(n = N, mean = mean(d$DCites)*p$df[RowSampled,"beta_d_cites_s"], sd = p$df[RowsSampled[x], "sigma"])
  ZPredictions2b <- rnorm(n = N, mean = mean(d$DCites)*p$df[RowSampled,"beta_d_cites_s"], sd = p$df[RowsSampled[x], "sigma"])
  
  Predictions <- XPredictions[[1]]@x + ZPredictions2a + ZPredictions2b
  
  PZero <- rbinom(n = N, size = 1, prob = logistic(Predictions))
  
  AllPredList[[x]] <- PZero

}

AllPredDF <- as.data.frame(AllPredList)
AllMammaldf$PredVirus <- apply(AllPredDF,1, function(a) a %>% mean)

AllMammaldf$PredVirusQ <- cut(AllMammaldf$PredVirus,
                                  breaks = c(-1:10/10),
                                  labels = c(0:10/10))

# Simulating the network #####

AllSimNets <- AllSimGraphs <- list()

for(i in 1:length(AllPredList)){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperMammals)] <- round(AllPredList[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                           union(AllMammaldf$Sp,AllMammaldf$Sp2))
  
  AllSimNets[[i]] <- as(AssMat, "dgCMatrix")
  
  AllSimGraphs[[i]] <- graph.adjacency(AssMat, mode = "undirected", diag = F)
  
  if(i%%10==0) print(i)
  
}

AllDegdf <- sapply(AllSimGraphs, function(a) degree(a)) #%>% as.data.frame
AllEigendf <- sapply(AllSimGraphs, function(a) eigen_centrality(a)$vector)# %>% as.data.frame

AllPredDegrees <- apply(AllDegdf, 1, mean)
AllPredEigen <- apply(AllEigendf, 1, mean)

Components <- sapply(AllSimGraphs, function(b) components(b)$no)

Cluster = sapply(AllSimGraphs, function(a) transitivity(a)) # all zero, don't bother

Betweenness1 = sapply(SimGraphs1, function(a) mean(betweenness(a)))
Betweenness1b = sapply(SimGraphs1b, function(a) mean(betweenness(a)))
Betweenness2 = sapply(SimGraphs2, function(a) mean(betweenness(a)))
Betweenness3 = sapply(SimGraphs3, function(a) mean(betweenness(a)))
Betweenness3b = sapply(SimGraphs3b, function(a) mean(betweenness(a)))

AllPrev = apply(AllPredDF, 2, Prev)
Prev1b = apply(PredDF1b, 2, Prev)
Prev2 = apply(PredDF2, 2, Prev)
Prev3 = apply(PredDF3, 2, Prev)
Prev3b = apply(PredDF3b, 2, Prev)

Closeness1 = sapply(SimGraphs1, function(a) mean(closeness(a)))
Closeness1b = sapply(SimGraphs1b, function(a) mean(closeness(a)))
Closeness2 = sapply(SimGraphs2, function(a) mean(closeness(a)))
Closeness3 = sapply(SimGraphs3, function(a) mean(closeness(a)))
Closeness3b = sapply(SimGraphs3b, function(a) mean(closeness(a)))