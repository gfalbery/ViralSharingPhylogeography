
# 2_Modelling mammal-mammal sharing ####

library(MCMCglmm); library(ggregplot); library(INLA)

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

# Modelling all mammal-mammal pairs ####

mf = 1

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),#,
                 G = list(G1 = list(V = diag(2), nu = 2)))

MC1 <- MCMCglmm(
  data = FinalHostMatrix,
  Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinCites + DomDom),
  rcov =~ idh(trait):units, 
  prior = prior.zi,
  random =~ us(trait):mm(Sp + Sp2),
  family = "zipoisson",
  pl = TRUE,
  nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin=3000*mf)

library(parallel)

setCores <- 6 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "FinalHostMatrix")
clusterExport(cl, "mf")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
ClusterMCStart <- Sys.time()
ZI_10runs <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(
    data = FinalHostMatrix,
    Virus ~ trait -1 + trait:(Space + Phylo + Space:Phylo + MinCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pl = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=3000*mf)
  
})

stopCluster(cl) # Stop running the parallel cluster
ClusterMCEnd <- Sys.time()

ClusterMCMC <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$Sol)))
ClusterMCMCv <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$VCV)))

SampleCluster <- as.mcmc(apply(bind_rows(ClusterMCMC), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv <- as.mcmc(apply(bind_rows(ClusterMCMCv), 2, function(a) as.mcmc(sample(a, 1000))))

ZIModel <- ZI_10runs[[1]]

ZIModel$Sol <- SampleCluster
ZIModel$VCV <- SampleClusterv

summary(ZIModel)