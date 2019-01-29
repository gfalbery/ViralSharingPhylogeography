
# Generating the network of viruses/hosts

load("ZI_runs.Rdata")

CountColumns <- list(1:7*2-1, 15:662)
ZIColumns <- list(1:7*2, 663:1310)

i = 1

N = nrow(FinalHostMatrix)

PredList1 <- list()

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

CountColumns <- list()

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, unlist(CountColumns)]
  ZIFXSample <- ClusterMCMC[RowSampled, unlist(ZIColumns)]
  
  XZMatrix <- cbind(ModelList[[i]]$X, ModelList[[i]]$Z)
  
  CountXZMatrix <- as.matrix(XZMatrix[1:N,unlist(CountColumns)])
  ZIXZMatrix <- as.matrix(XZMatrix[(N+1):(2*N),unlist(ZIColumns)])
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput)
  PCount <- exp(CountOutput)*(1-PZero)
  
  PredList1[[x]] <- PCount
  
}

PredDF1 <- as.data.frame(PredList1)
MeanPredictions <- apply(PredDF1,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus1 <- MeanPredictions

ggplot(FinalHostMatrix, aes(Virus, PredVirus1)) + coord_fixed() + geom_point() + geom_smooth(method = lm)

# Trying it without random effects ####

i = 2

PredList2 <- list()

for(x in 1:1000){
  
  print(x)
  
  RowSampled <- sample(1:nrow(ModelList[[i]]$Sol), 1)
  
  CountFXSample <- ModelList[[i]]$Sol[RowSampled, CountColumns[[1]]]
  ZIFXSample <- ModelList[[i]]$Sol[RowSampled, ZIColumns[[1]]]
  
  XZMatrix <- cbind(ModelList[[i]]$X, ModelList[[i]]$Z)
  
  CountXZMatrix <- as.matrix(XZMatrix[,CountColumns[[1]]])
  ZIXZMatrix <- as.matrix(XZMatrix[,ZIColumns[[1]]])
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  CountOutput <- CountOutput[1:N]
  ZIOutput <- ZIOutput[(N+1):(2*N)]
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput)
  PCount <- exp(CountOutput)*(1-PZero)
  
  PredList2[[x]] <- PCount
  
}

PredDF2 <- as.data.frame(PredList2)
ModePredictions <- apply(PredDF2,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus2 <- ModePredictions

ggplot(FinalHostMatrix, aes(PredVirus, Virus)) + coord_fixed() + geom_point() + geom_smooth(method = lm)
ggplot(FinalHostMatrix, aes(PredVirus, PredVirus2)) + coord_fixed() + geom_point() + geom_smooth(method = lm)

# Simulating the networks #####

SimNets1 <- SimGraphs1 <- vector("list", length(PredList1))

for(i in 1:length(PredList1)){
  
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

Degdf <- sapply(SimGraphs1, function(a) degree(a)) %>% as.data.frame

Degdflong <- reshape2::melt(t(Degdf)) %>% rename(Sp = Var2, Degree = value)

a = 850

ggplot(Degdflong[1:a,], aes(Sp, Degree)) + geom_violin(aes(colour = Sp)) +
  geom_point(data = Hosts[Hosts$Sp%in%Degdflong[1:a,"Sp"],], aes(Sp, Degree)) + 
  facet_wrap(~Sp, scales = "free")

PredDegrees <- apply(Degdf, 1, mean)

Hosts$PredDegree1 <- PredDegrees[as.character(Hosts$Sp)]

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

Degdf2 <- lapply(SimGraphs2[1:609], function(a) degree(a)) %>% as.data.frame

Degdflong2 <- reshape2::melt(t(Degdf)) %>% rename(Sp = Var2, Degree = value)

a = 850

ggplot(Degdflong2[1:a,], aes(Sp, Degree)) + geom_violin(aes(colour = Sp)) +
  geom_point(data = Hosts[Hosts$Sp%in%Degdflong2[1:a,"Sp"],], aes(Sp, Degree)) + 
  facet_wrap(~Sp, scales = "free")

PredDegrees <- apply(Degdf, 1, mean)

Hosts$PredDegree2 <- PredDegrees[as.character(Hosts$Sp)]

ggplot(Hosts, aes(Degree, PredDegree2)) + geom_point() + geom_smooth()

