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

ZI_runs2 <- parallel::mclapply(1:40, function(i) {
  if(i <= 10) {
    
    return(MCMCglmm( # Running full matrix with simple random effect of row-species
      data = FinalHostMatrix2,
      Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi,
      random =~ us(trait):Sp,
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf)
    )
    
  } else if (i > 10 & i <= 20) {
    
    return(MCMCglmm( # Running full model with spacequantile effect
      data = FinalHostMatrix,
      Virus ~ trait -1 + trait:(Space + Phylo2 + SpaceQuantile:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi,
      random =~ us(trait):mm(Sp + Sp2),
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf)
    )
    
  } else if (i > 20 & i <= 30) {
    
    return(MCMCglmm( # Only spatial overlap with quantiles
      data = FinalHostMatrixNoSpace,
      Virus ~ trait -1 + trait:(Space + Phylo2 + SpaceQuantile:Phylo2 + MinDCites + DomDom),
      rcov =~ idh(trait):units, 
      prior = prior.zi,
      random =~ us(trait):mm(Sp + Sp2),
      family = "zipoisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf))
    
  } else if (i > 30 & i <= 40) {
    
    save(MCMCglmm( # Only spatial overlap with quantiles
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

# Save file

save(ZI_runs2, file = "ZI_runs2.Rdata")
save.image()
