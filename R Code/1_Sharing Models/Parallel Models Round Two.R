# Parallel models 2: running some less problematic stuff ####

# Run source code

rm(list = ls())

source("R Code/00_Master Code.R")

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel); library(dplyr)

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostThemselves <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo2
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrixNoSpace <- FinalHostMatrix %>% filter(Space>0)

FinalHostMatrixNoSpace <- droplevels(FinalHostMatrix[FinalHostMatrix$Space>0,])
FinalHostMatrixNoSpace$Sp <- factor(FinalHostMatrixNoSpace$Sp, levels = union(FinalHostMatrixNoSpace$Sp,FinalHostMatrixNoSpace$Sp2))
FinalHostMatrixNoSpace$Sp2 <- factor(FinalHostMatrixNoSpace$Sp2, levels = union(FinalHostMatrixNoSpace$Sp,FinalHostMatrixNoSpace$Sp2))

FinalHostMatrix2 <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix2$Phylo <- FinalHostMatrix2$Phylo2
FinalHostMatrix2$MinDCites <- log(FinalHostMatrix2$MinDCites + 1)
FinalHostMatrix2NoSpace <- FinalHostMatrix2 %>% filter(Space>0)

gp <- gelman.prior(Virus ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom, data = FinalHostMatrix2)

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),
                 G = list(G1 = list(V = diag(2), nu = 2)),
                 B = list(V = diag(14)*10^8, mu = rep(0,14)))

prior.zi$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp

prior.zi2 <- list(R = list(V = diag(2), nu = 0, fix = 2),
                  #G = list(G1 = list(V = diag(2), nu = 2)),
                  B = list(V = diag(14)*10^8, mu = rep(0,14)))

prior.zi2$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp

# Modelling all mammal-mammal pairs ####

mf = 15


prior.pois <- list(R = list(V = diag(1), nu = 0),
                 G = list(G1 = list(V = diag(1), nu = 2)))

prior.pois2 <- list(R = list(V = diag(1), nu = 0))

# Trying a Poisson model ####

parallel::mclapply(1:40, function(i) {
  if(i <= 10) {
    
    saveRDS(MCMCglmm(
      data = FinalHostMatrix,
      Virus ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      rcov =~ idh(trait):units, 
      prior = prior.pois,
      random =~ mm(Sp + Sp2),
      family = "poisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Poisson Model ",i, ".Rdata"))
    
  } else if (i > 10 & i <= 20) {
    
    saveRDS(MCMCglmm(
      data = FinalHostMatrix,
      Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.pois2,
      #random =~ us(trait):mm(Sp + Sp2),
      family = "poisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Poisson Model ",i, ".Rdata"))
  }
})

parallel::mclapply(1:40, function(i) {
  if(i <= 10) {
    
    saveRDS(MCMCglmm( # Running full matrix with simple random effect of row-species
      data = FinalHostMatrix2,
      Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi,
      random =~ us(trait):Sp,
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Model",i, ".Rdata"))
    
  } else if (i > 10 & i <= 20) {
    
    saveRDS(MCMCglmm( # Running full model with spacequantile effect
      data = FinalHostMatrix,
      Virus ~ trait -1 + trait:(Space + Phylo2 + SpaceQuantile:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi,
      random =~ us(trait):mm(Sp + Sp2),
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Model",i, ".Rdata"))
    
  } else if (i > 20 & i <= 30) {
    
    saveRDS(MCMCglmm( # Only spatial overlap with quantiles
      data = FinalHostMatrixNoSpace,
      Virus ~ trait -1 + trait:(Space + Phylo2 + SpaceQuantile:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi,
      random =~ us(trait):mm(Sp + Sp2),
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Model",i, ".Rdata"))
    
  } else if (i > 30 & i <= 40) {
    
    saveRDS(MCMCglmm( # Only spatial overlap with quantiles
      data = FinalHostMatrixNoSpace,
      Virus ~ trait -1 + trait:(Space + Phylo2 + SpaceQuantile:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi2,
      #random =~ us(trait):mm(Sp + Sp2),
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Model",i, ".Rdata"))
  }}, mc.cores = 40)

# Trying Binomial test ####

FinalHostMatrix$VirBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

prior.zi3 <- list(R = list(V = diag(1), nu = 0, fix = 1),
                  G = list(G1 = list(V = diag(1), nu = 2)))

prior.zi4<- list(R = list(V = diag(1), nu = 0, fix = 1))

parallel::mclapply(1:40, function(i) {
  if(i <= 10) {
    
    saveRDS(MCMCglmm( # Running binomial model with simple random effect of row-species
      data = FinalHostMatrix,
      VirBinary ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi3,
      random =~ us(trait):mm(Sp+Sp2),
      family = "categorial",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, 
      trunc = T), file = paste0("BinModel",i, ".Rdata"))
    
  } else if (i > 10 & i <= 20) {
    
    saveRDS(MCMCglmm( # Running full model with spacequantile effect
      data = FinalHostMatrix,
      VirBinary ~ Space + Phylo2 + SpaceQuantile:Phylo2 + MinDCites + DomDom,
      prior = prior.zi4,
      family = "categorical",
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, 
      trunc = T), file = paste0("BinModel",i, ".Rdata"))}
  
}) 