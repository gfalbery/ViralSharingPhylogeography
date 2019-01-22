
# 2_Modelling mammal-mammal sharing ####

library(MCMCglmm); library(ggregplot); library(INLA)

# Modelling all mammal-mammal pairs ####

mf = 5

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),#,
                 G = list(G1 = list(V = diag(2), nu = 2)))

library(parallel)

setCores <- 6 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "HostMatrixdf")
clusterExport(cl, "UpperHosts")
clusterExport(cl, "mf")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
ClusterMCStart <- Sys.time()
ZI_10runs <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(#data = HostMatrixdf[-HostThemselves,], 
    data = HostMatrixdf[-HostThemselves,],
    Virus ~ trait -1 + trait:(Space + Phylo + SpaceQuantile:Phylo + Cites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):Sp,
    family = "zipoisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)
  
})

stopCluster(cl) # Stop running the parallel cluster
ClusterMCEnd <- Sys.time()

ClusterMCMC <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$Sol)))
ClusterMCMCv <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$VCV)))

SampleCluster <- as.mcmc(apply(bind_rows(ClusterMCMC), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv <- as.mcmc(apply(bind_rows(ClusterMCMCv), 2, function(a) as.mcmc(sample(a, 1000))))

SampleCluster %>% 
  apply(2, function(a) {
    qplot(1:1000, a, geom = "line") + labs(x = "Iteration", y = "Effect", title = colnames(a))
  }) %>% 
  arrange_ggplot2

ZIModel <- ZI_10runs[[1]]

ZIModel$Sol <- SampleCluster
ZIModel$VCV <- SampleClusterv
summary(ZIModel)

# Modelling only mammal-mammal pairs that share space ####

mf = 5

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),#,
                 G = list(G1 = list(V = diag(2), nu = 2)))

library(parallel)

OverlapHostMatrixdf <- droplevels(HostMatrixdf[-HostThemselves&HostMatrixdf$Space>0,])

setCores <- 6 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "OverlapHostMatrixdf")
clusterExport(cl, "mf")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
ClusterMCStart2 <- Sys.time()
ZI_10runs2 <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(#data = HostMatrixdf[-HostThemselves,], 
    data = OverlapHostMatrixdf,
    Virus ~ trait -1 + trait:(Space + Phylo + SpaceQuantile:Phylo + Cites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):Sp,
    family = "zipoisson",
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)
  
})

stopCluster(cl) # Stop running the parallel cluster
ClusterMCEnd2 <- Sys.time()

ClusterMCMC2 <- ZI_10runs2 %>% lapply(function(a) as.data.frame(as.matrix(a$Sol)))
ClusterMCMCv2 <- ZI_10runs2 %>% lapply(function(a) as.data.frame(as.matrix(a$VCV)))

SampleCluster2 <- as.mcmc(apply(bind_rows(ClusterMCMC2), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv2 <- as.mcmc(apply(bind_rows(ClusterMCMCv2), 2, function(a) as.mcmc(sample(a, 1000))))

SampleCluster %>% 
  apply(2, function(a) {
    qplot(1:1000, a, geom = "line") + labs(x = "Iteration", y = "Effect", title = colnames(a))
  }) %>% 
  arrange_ggplot2

ZIModelNoSpace <- ZI_10runs2[[1]]

ZIModelNoSpace$Sol <- SampleCluster2
ZIModelNoSpace$VCV <- SampleClusterv2
summary(ZIModelNoSpace)
