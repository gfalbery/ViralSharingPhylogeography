
i = 3

mc <- BinModelList[1:10 + 10*(i-1)] %>% lapply(function(a) a$Sol[,1:7])
#mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$VCV)
mc <- do.call(mcmc.list, mc)
#par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
#gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(7,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

N = nrow(FinalHostMatrix)

# Going from the full model ####

PredList4 <- list()

ClusterMCMC <- BinModelList[1:10 + 10*(i-1)] %>% 
  lapply(function(a) as.data.frame(a$Sol)) %>% 
  bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- cbind(BinModelList[[i*10]]$X, BinModelList[[i*10]]$Z) %>% 
  as.matrix %>% as("dgCMatrix")

Columns <- list(1:ncol(BinModelList[[i*10]]$X),(ncol(BinModelList[[i*10]]$X)+1):ncol(XZMatrix))

for(x in 1:1000){
  if(x%%10==0) print(x)
  RowSampled <- RowsSampled[x]
  FXSample <- ClusterMCMC[RowSampled, unlist(Columns)]
  FXSample <- FXSample# + 0.04
  Output <- c(FXSample %*% t(XZMatrix))
  ProbVector <- Output[[1]]@x
  PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector))
  #PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector)*Prev(FinalHostMatrix$VirusBinary)/mean(invlogit(ProbVector)))
  PredList4[[x]] <- PZero
}

PredDF4 <- as.data.frame(PredList4)
names(PredDF4) <- paste("Rep",1:1000)
FinalHostMatrix$PredVirus4 <- apply(PredDF1,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus4Q <- cut(FinalHostMatrix$PredVirus4,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Simulating Networks ####

SimNets4 <- SimGraphs4 <- list()

for(i in 1:length(PredList4)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList4[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1 <- sapply(SimGraphs1, function(a) degree(a)) %>% as.data.frame
Eigendf1 <- sapply(SimGraphs1, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1 <- apply(Degdf1, 1, mean)
PredEigen1 <- apply(Eigendf1, 1, mean)

Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]
Hosts$PredEigen1 <- PredEigen1[as.character(Hosts$Sp)]

# Without random effects, same model ####

i = 3

PredList4b <- list()

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- BinModelList[[i*10]]$X %>% 
  as.matrix %>% as("dgCMatrix")

for(x in 1:1000){
  if(x%%10==0) print(x)
  RowSampled <- RowsSampled[x]
  FXSample <- ClusterMCMC[RowSampled, Columns[[1]]]
  Output <- c(FXSample %*% t(XZMatrix))
  ProbVector <- Output[[1]]@x
  PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector))
  #PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector)*Prev(FinalHostMatrix$VirusBinary)/mean(invlogit(ProbVector)))
  PredList4b[[x]] <- PZero
}

PredDF4b <- as.data.frame(PredList4b)
FinalHostMatrix$PredVirus1b <- apply(PredDF1b[,1:1000], 1, function(a) a %>% mean)

# Simulating the networks #####

SimNets1b <- SimGraphs1b <- list()

for(i in 1:length(PredList4b)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList4b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1b[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1b[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1b <- sapply(SimGraphs1b, function(a) degree(a)) %>% as.data.frame
Eigendf1b <- sapply(SimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1b <- apply(Degdf1b, 1, mean)
PredEigen1b <- apply(Eigendf1b, 1, mean)

Hosts$PredDegree1b <- PredDegrees1b[as.character(Hosts$Sp)]
Hosts$PredEigen1b <- PredEigen1b[as.character(Hosts$Sp)]

FinalHostMatrix$PredVirus1bQ <- cut(FinalHostMatrix$PredVirus1b,
                                    breaks = c(-1:10/10),
                                    labels = c(0:10/10))
