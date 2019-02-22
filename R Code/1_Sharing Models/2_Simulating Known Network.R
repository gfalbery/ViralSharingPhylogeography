
# STAN Model Output ####

library(rstan); library(reskew); library(ggregplot); library(parallel)

# Doing the simulating #####

PredList1 <- list()

Predictions <- predict(BAMList[[1]], newdata = DataList[[1]])

N = nrow(FinalHostMatrix)

PredList1 <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
  
  BinPred <- rbinom(n = N,
                    prob = logistic(Predictions),
                    size  = 1)
  
  BinPred
  
}, mc.cores = 10)

PredDF1 <- data.frame(PredList1)
FinalHostMatrix$PredVirus1 <- apply(PredDF1, 1, mean)
FinalHostMatrix$PredVirus1Q <- cut(FinalHostMatrix$PredVirus1,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Simulating without the random effect ####

SimNets1 <- mclapply(1:length(PredList1), function(i){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[lower.tri(AssMat)] <- round(PredList1[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

SimGraphs1 <- mclapply(1:length(PredList1), function(i){
  
  graph.adjacency(SimNets1[[i]], mode = "undirected")
  
}, mc.cores = 10)

Degdf1 <- sapply(SimGraphs1, function(a) degree(a)) %>% as.data.frame
PredDegrees1 <- apply(Degdf1, 1, mean)
Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]

Eigendf1 <- sapply(SimGraphs1, function(a) eigen_centrality(a)$vector) %>% as.data.frame
PredEigen1 <- apply(Eigendf1, 1, mean)
Hosts$PredEigen1 <- PredEigen1[as.character(Hosts$Sp)]

# Getting network-level stats

SimGraphList <- list(SimGraphs1, SimGraphs1)# , SimGraphs2, SimGraphs3, SimGraphs3b)

Components <- lapply(SimGraphList, function(a) sapply(a, function(b) components(b)$no))
ComponentSizes <- lapply(SimGraphList, function(a) sapply(a, function(b) components(b)$no))
lapply(SimGraphs1[which(Components[[2]]==2)], function(a) which(components(a)$membership==2))

Degrees <- lapply(SimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(SimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1 = sapply(SimGraphs1, function(a) transitivity(a))

#Betweenness1 = sapply(SimGraphs1, function(a) betweenness(a))
#Betweenness1 = sapply(SimGraphs1, function(a) betweenness(a))

Prev1 = sapply(PredList1, Prev)

#Closeness1 = sapply(SimGraphs1, function(a) closeness(a))
#Closeness1b = sapply(SimGraphs1b, function(a) closeness(a))

