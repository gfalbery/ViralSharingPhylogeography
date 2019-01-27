
# Summarising parallel models ####

mc <- ZI_runs[1:10] %>% lapply(function(a) a$Sol[,1:14])
mc <- do.call(mcmc.list, mc)
par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(14,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

mc <- ZI_runs[1:10+10] %>% lapply(function(a) a$Sol[,1:14])
mc <- do.call(mcmc.list, mc)
par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(14,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

ClusterMCMC <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$Sol[,1:14]))) %>% bind_rows
ClusterMCMCv <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$VCV))) %:% bind_rows

SampleCluster <- as.mcmc(apply(bind_rows(ClusterMCMC), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv <- as.mcmc(apply(bind_rows(ClusterMCMCv), 2, function(a) as.mcmc(sample(a, 1000))))


FullZIModel <- ZI_runs[[1]]
ZINoGModel <- ZI_runs[[11]]
OverlapZIModel <- ZI_runs[[21]]
OverlapNoGZIModel <- ZI_runs[[31]]

ModelList <- list(FullZIModel, ZINoGModel, OverlapZIModel, OverlapNoGZIModel)

for(i in 1:4){
  
  ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol[,1:14]))) %>% bind_rows
  ClusterMCMCv <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$VCV))) %>% bind_rows
  
  rows <- sample(1:dim(ClusterMCMC)[1], 1000)
  
  SampledSol <- ClusterMCMC[rows,]
  SampledVCV <- ClusterMCMCv[rows,]
  
  ModelList[[i]]$Sol <- SampledSol
  ModelList[[i]]$VCV <- SampledVCV
  
}




