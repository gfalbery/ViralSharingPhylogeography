
# Generating the network of viruses/hosts

library(MCMCglmm)

load("ZI_runs.Rdata")

CountColumns <- list(1:7*2-1, 15:662)
ZIColumns <- list(1:7*2, 663:1310)

i = 1

N = nrow(FinalHostMatrix)

PredList1 <- list()

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

logit <- function(a) exp(a)/(1 + exp(a))

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, unlist(CountColumns)]
  ZIFXSample <- ClusterMCMC[RowSampled, unlist(ZIColumns)]
  
  XZMatrix <- cbind(ZI_runs[[(i-1)*10+1]]$X, ZI_runs[[(i-1)*10+1]]$Z)
  
  CountXZMatrix <- XZMatrix[1:N,unlist(CountColumns)] #%>% as.matrix 
  ZIXZMatrix <- XZMatrix[(N+1):(2*N),unlist(ZIColumns)] #%>% as.matrix
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput[[1]]@x)
  PCount <- exp(CountOutput[[1]]@x)*(1-PZero)
  
  PredList1[[x]] <- PCount
  
}

PredDF1 <- as.data.frame(PredList1)
MeanPredictions <- apply(PredDF1,1, function(a) a %>% mean)
FinalHostMatrix$PredVirus1 <- MeanPredictions

# Without random effects, same model ####

PredList1b <- list()

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, CountColumns[[1]]]
  ZIFXSample <- ClusterMCMC[RowSampled, ZIColumns[[1]]]
  
  XZMatrix <- cbind(ZI_runs[[(i-1)*10+1]]$X, ZI_runs[[(i-1)*10+1]]$Z)
  
  CountXZMatrix <- XZMatrix[1:N,CountColumns[[1]]] #%>% as.matrix
  ZIXZMatrix <- XZMatrix[(N+1):(2*N),ZIColumns[[1]]] #%>% as.matrix
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput[[1]]@x)
  PCount <- exp(CountOutput[[1]]@x)*(1-PZero)
  
  PredList1b[[x]] <- PCount
  
}

# Trying it without random effects ####

i = 2

PredList2 <- list()

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, CountColumns[[1]]]
  ZIFXSample <- ClusterMCMC[RowSampled, ZIColumns[[1]]]
  
  XZMatrix <- cbind(ZI_runs[[(i-1)*10+1]]$X, ZI_runs[[(i-1)*10+1]]$Z)
  
  CountXZMatrix <- XZMatrix[1:N,CountColumns[[1]]] #%>% as.matrix
  ZIXZMatrix <- XZMatrix[(N+1):(2*N),ZIColumns[[1]]] #%>% as.matrix
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput[[1]]@x)
  PCount <- exp(CountOutput[[1]]@x)*(1-PZero)
  
  PredList2[[x]] <- PCount
  
}

PredDF2 <- as.data.frame(PredList2)
ModePredictions <- apply(PredDF2,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus2 <- ModePredictions

# Simulating the networks #####

SimNets1 <- SimGraphs1 <- vector("list", length(PredList1))

for(i in 1:length(PredList1)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1[[i]] <- AssMat
  
  SimGraphs1[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1 <- sapply(SimGraphs1, function(a) degree(a)) %>% as.data.frame
Eigendf1 <- sapply(SimGraphs1, function(a) eigen_centrality(a)$vector) %>% as.data.frame

Degdflong1 <- reshape2::melt(t(Degdf1)) %>% rename(Sp = Var2, Degree = value)

a = 850

ggplot(Degdflong1[1:a,], aes(Sp, Degree)) + geom_violin(aes(colour = Sp)) +
  geom_point(data = Hosts[Hosts$Sp%in%Degdflong1[1:a,"Sp"],], aes(Sp, Degree)) + 
  facet_wrap(~Sp, scales = "free")

PredDegrees1 <- apply(Degdf1, 1, mean)

Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]

# Simulating the networks #####

SimNets1b <- SimGraphs1b <- vector("list", length(PredList1b))

for(i in 1:length(PredList1b)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1b[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1b[[i]] <- AssMat
  
  SimGraphs1b[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1b <- sapply(SimGraphs1b, function(a) degree(a)) %>% as.data.frame
Eigendf1b <- sapply(SimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

Degdflong1b <- reshape2::melt(t(Degdf1b)) %>% rename(Sp = Var2, Degree = value)

a = 850

ggplot(Degdflong1b[1:a,], aes(Sp, Degree)) + geom_violin(aes(colour = Sp)) +
  geom_point(data = Hosts[Hosts$Sp%in%Degdflong1b[1:a,"Sp"],], aes(Sp, Degree)) + 
  facet_wrap(~Sp, scales = "free")

PredDegrees1b <- apply(Degdf1b, 1, mean)

Hosts$PredDegree1b <- PredDegrees1b[as.character(Hosts$Sp)]

# Trying sans random effect ####

SimNets2 <- SimGraphs2 <- vector("list", length(PredList2))

for(i in 1:length(PredList2)){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList2[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets2[[i]] <- AssMat
  
  SimGraphs2[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
  if(i%%10==0) print(i)
  
}

Degdf2 <- sapply(SimGraphs2, function(a) degree(a)) %>% as.data.frame
Eigendf2 <- sapply(SimGraphs2[1:729], function(a) eigen_centrality(a)$vector) %>% as.data.frame

Degdflong2 <- reshape2::melt(t(Degdf2)) %>% rename(Sp = Var2, Degree = value)

a = 850

ggplot(Degdflong2[1:a,], aes(Sp, Degree)) + geom_violin(aes(colour = Sp)) +
  geom_point(data = Hosts[Hosts$Sp%in%Degdflong2[1:a,"Sp"],], aes(Sp, Degree)) + 
  facet_wrap(~Sp, scales = "free")

PredDegrees2 <- apply(Degdf2, 1, mean)

Hosts$PredDegree2 <- PredDegrees2[as.character(Hosts$Sp)]

ggplot(Hosts, aes(PredDegree1b, PredDegree2)) + geom_point() + geom_smooth()

GGally::ggpairs(Hosts[,c("Degree","PredDegree1","PredDegree1b","PredDegree2")],
                lower = list(continuous = "smooth"))

apply(Hosts[,c("Degree","PredDegree1","PredDegree1b","PredDegree2")],2,function(a) table(a>0))

PredEigen1 <- apply(Eigendf1,1,mean)
PredEigen1b <- apply(Eigendf1b,1,mean)
PredEigen2 <- apply(Eigendf2,1,mean)

Hosts[,c("Eigen1","Eigen1b","Eigen2")] <- cbind(PredEigen1[as.character(Hosts$Sp)],
                                                PredEigen1b[as.character(Hosts$Sp)],
                                                PredEigen2[as.character(Hosts$Sp)])

GGally::ggpairs(Hosts[,c("Eigenvector","Eigen1","Eigen1b","Eigen2")],
                lower = list(continuous = "smooth"))


# Comparing differences ####

BeforeHPD <- cbind(apply(Degdf1[1:(nrow(Degdf1)/2),],1, function(a) HPDinterval(as.mcmc(a))[1]),
                   apply(Degdf1[1:(nrow(Degdf1)/2),],1, function(a) HPDinterval(as.mcmc(a))[2])) %>% as.data.frame

AfterHPD <- cbind(apply(Degdf1b[1:(nrow(Degdf1b)/2),],1, function(a) HPDinterval(as.mcmc(a))[1]),
                  apply(Degdf1b[1:(nrow(Degdf1b)/2),],1, function(a) HPDinterval(as.mcmc(a))[2])) %>% as.data.frame
BeforeHPD$When <- "Before"
AfterHPD$When <- "After"
BeforeHPD$Sp <- FHN
AfterHPD$Sp <- FHN
HPDComp <- rbind(BeforeHPD, AfterHPD)
HPDComp2 <- cbind(BeforeHPD, AfterHPD)
names(HPDComp2) <- paste(names(HPDComp2),rep(1:2,each = 4), sep = ".")







