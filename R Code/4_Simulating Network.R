
# Generating the network of viruses/hosts

load("ZI_runs.Rdata")

FullZIModel <- ZI_runs[[1]]
ZINoGModel <- ZI_runs[[11]]
OverlapZIModel <- ZI_runs[[21]]
OverlapNoGZIModel <- ZI_runs[[31]]

ModelList <- list(FullZIModel, ZINoGModel, OverlapZIModel, OverlapNoGZIModel)

i = 1

PredList1 <- list()

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

for(x in 1:1000){
  
  RowSampled <- sample(1:nrow(ClusterMCMC), 1)
  
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

PredDF1 <- as.data.frame(PredList)
ModePredictions <- apply(PredDF1,1, function(a) a %>% as.mcmc %>% posterior.mode)
MeanPredictions <- apply(PredDF1,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus <- ModePredictions

ggplot(FinalHostMatrix, aes(PredVirus, Virus)) + coord_fixed() + geom_point() + geom_smooth(method = lm)

# Trying it without random effects ####

i = 2

PredList2 <- list()

for(x in 1:1000){
  
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
