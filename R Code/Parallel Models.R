
# Run source code

rm(list = ls())

source("R Code/00_Master Code.R")

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel)

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo2
FinalHostMatrixNoSpace <- FinalHostMatrix %>% filter(Space>0)

gp <- gelman.prior(Virus ~ Space + Phylo2 + Space:Phylo2 + MinCites + DomDom, data = FinalHostMatrix)

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),
                 G = list(G1 = list(V = diag(2), nu = 2)),
                 B = list(V = diag(16)*10^8, mu = rep(0,16)))

prior.zi$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp

prior.zi2 <- list(R = list(V = diag(2), nu = 0, fix = 2),
                  #G = list(G1 = list(V = diag(2), nu = 2)),
                  B = list(V = diag(16)*10^8, mu = rep(0,16)))

prior.zi2$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp

# Modelling all mammal-mammal pairs ####

mf = 15

setCores <- 20 # use detectCores() by itself if you want all CPUs

cl <- makeCluster(getOption("cl.cores", setCores))

cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) # load the MCMCglmm package within the cluster

clusterExport(cl, "prior.zi") # Import each object that's necessary to run the function
clusterExport(cl, "prior.zi2") # Import each object that's necessary to run the function
clusterExport(cl, "FinalHostMatrix")
clusterExport(cl, "mf")

ZI_runs <- parLapply(cl = cl, 1:10, function(i) {
  
  ZI_runsFull <- 
    
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
  
  # All mammals no G ####
  
  ZI_runsNoG <- 
    
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
  
  # Only overlapping mammals ####
  
  ZI_runsNoSpace <- 
    
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
  
  # Only overlapping mammals no G ####
  
  ZI_runsNoSpaceNoG <- 
    
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
  
  return(list(ZI_runsFull, ZI_runsNoG, ZI_runsNoSpace, ZI_runsNoSpaceNoG))
  
})

stopCluster(cl) # Stop running the parallel cluster

# Save file

save(ZI_runs, file = "ZI_runs.Rdata")
save.image()
