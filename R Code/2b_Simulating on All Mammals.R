
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

AllPredDF <- as.data.frame(AllPredList) %>% as.matrix %>% unname %>% as("dgCMatrix")
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

Sgs <- parallel::mclapply(1:10, function(i){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperMammals)] <- round(AllPredList[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                           union(AllMammaldf$Sp,AllMammaldf$Sp2))
  
  list(
  
  as(AssMat, "dgCMatrix"),
  
  graph.adjacency(AssMat, mode = "undirected", diag = F)
  
  ) %>% return
  #if(i%%10==0) print(i)
  
}, mc.cores = 10)

AllSims <- map(Sgs, function(a) a[[1]])
AllSimgs <- map(Sgs, function(a) a[[2]])

AllPrev = apply(AllPredDF, 2, Prev)

AllDegdf <- sapply(AllSimGraphs, function(a) degree(a)) #%>% as.data.frame
AllEigendf <- sapply(AllSimGraphs, function(a) eigen_centrality(a)$vector)# %>% as.data.frame

AllPredDegrees <- apply(AllDegdf, 1, mean)
AllPredEigen <- apply(AllEigendf, 1, mean)

Components <- sapply(AllSimGraphs, function(b) components(b)$no)
Cluster = sapply(AllSimGraphs, function(a) transitivity(a)) # all zero, don't bother

Betweenness1 = sapply(AllSimGraphs, function(a) mean(betweenness(a)))

AllCloseness = sapply(AllSimGraphs, function(a) mean(closeness(a)))


Hosts[,"AllPredDegree"] <- AllPredDegrees[as.character(Hosts$Sp)]

