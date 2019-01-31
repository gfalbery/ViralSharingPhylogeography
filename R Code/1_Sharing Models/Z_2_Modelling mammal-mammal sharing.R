
# 2_Modelling mammal-mammal sharing ####

library(MCMCglmm); library(ggregplot); library(INLA)

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo2
FinalHostMatrixNoSpace <- FinalHostMatrix %>% filter(Space>0)

# Modelling all mammal-mammal pairs ####

library(parallel)

mf = 6

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
    Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf)
  
})

stopCluster(cl) # Stop running the parallel cluster
ClusterMCEnd <- Sys.time()

mc <- ZI_10runs %>% lapply(function(a) a$Sol[,1:14])
mc <- do.call(mcmc.list, mc)
par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(14,2), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

ClusterMCMC <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$Sol[,1:14])))
ClusterMCMCv <- ZI_10runs %>% lapply(function(a) as.data.frame(as.matrix(a$VCV)))

SampleCluster <- as.mcmc(apply(bind_rows(ClusterMCMC), 2, function(a) as.mcmc(sample(a, 1000))))
SampleClusterv <- as.mcmc(apply(bind_rows(ClusterMCMCv), 2, function(a) as.mcmc(sample(a, 1000))))

ZIModel <- ZI_10runs[[1]]

ZIModel$Sol <- SampleCluster
ZIModel$VCV <- SampleClusterv

summary(ZIModel)

ZISols <- summary(ZIModel)$solutions %>% as.data.frame %>% mutate(
  Component = rep(c("Count", "ZI"), dim(ZISols)[1]/2),
  Variable = rep(c("Intercept", "Space", "Phylogeny", "Citations", "DomDom", "DomWild", 
                   "Phylo:Space"),
                 each = 2),
  Name = paste(Component, Variable, sep = ":")
) %>% rename(Lower = "l-95% CI", Upper = "u-95% CI", Estimate = "post.mean")

ZISols$Estimate[ZISols$Component=="ZI"] <- -ZISols$Estimate[ZISols$Component=="ZI"]
ZISols[ZISols$Component=="ZI", c("Lower", "Upper")] <- -ZISols[ZISols$Component=="ZI", c("Upper","Lower")]

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = Component)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + coord_flip() + 
  facet_wrap(~Component, scale = "free") +
  ggtitle("Zero-Inflated Model Output") 

# Parallelising all model runs #####

gp <- gelman.prior(Virus ~ Space + Phylo2 + Space:Phylo2 + MinCites + DomDom, data = FinalHostMatrix)

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),
                  G = list(G1 = list(V = diag(2), nu = 2)),
                  B = list(V = diag(16)*10^8, mu = rep(0,16)))

prior.zi$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp

prior.zi2 <- list(R = list(V = diag(2), nu = 0, fix = 2),
                  #G = list(G1 = list(V = diag(2), nu = 2)),
                  B = list(V = diag(16)*10^8, mu = rep(0,16)))

prior.zi2$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp


# make the cluster

mf = 10

library(parallel)

setCores <- 20 # use detectCores() by itself if you want all CPUs

cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "prior.zi2") # Import each object that's necessary to run the function
clusterExport(cl, "FinalHostMatrix")
clusterExport(cl, "mf")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
Time0 <- Sys.time()

# Full model ####

ZI_runsFull <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(
    data = FinalHostMatrix,
    Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf)
  
})

Time1 <- Sys.time()

# All mammals no G ####

ZI_runsNoG <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(
    data = FinalHostMatrix,
    Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi2,
    #random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf)
  
})

Time2 <- Sys.time()

# Only overlapping mammals ####

ZI_runsNoSpace <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(
    data = FinalHostMatrix %>% filter(Space>0),
    Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf)
  
})

Time3 <- Sys.time()

# Only overlapping mammals no G ####

ZI_runsNoG <- parLapply(cl = cl, 1:10, function(i) {
  
  MCMCglmm(
    data = FinalHostMatrix %>% filter(Space>0),
    Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi2,
    #random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf)
  
})

stopCluster(cl) # Stop running the parallel cluster

Time4 <- Sys.time()
