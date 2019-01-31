# Parallel Poisson models ####


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
FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))

FinalHostMatrixNoSpace <- droplevels(FinalHostMatrix[FinalHostMatrix$Space>0,])
FinalHostMatrixNoSpace$Sp <- factor(FinalHostMatrixNoSpace$Sp, levels = union(FinalHostMatrixNoSpace$Sp,FinalHostMatrixNoSpace$Sp2))
FinalHostMatrixNoSpace$Sp2 <- factor(FinalHostMatrixNoSpace$Sp2, levels = union(FinalHostMatrixNoSpace$Sp,FinalHostMatrixNoSpace$Sp2))

# Modelling all mammal-mammal pairs ####

mf = 15

prior.pois <- list(R = list(V = diag(1), nu = 0),
                   G = list(G1 = list(V = diag(1), nu = 2)))

prior.pois2 <- list(R = list(V = diag(1), nu = 0))

# Trying a Poisson model ####

parallel::mclapply(1:20, function(i) {
  if(i <= 10) {
    
    saveRDS(MCMCglmm(
      data = FinalHostMatrix,
      Virus ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      prior = prior.pois,
      random =~ mm(Sp + Sp2),
      family = "poisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Model Runs/Poisson Model ",i, ".Rdata"))
    
  } else if (i > 10 & i <= 20) {
    
    saveRDS(MCMCglmm(
      data = FinalHostMatrix,
      Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom),
      prior = prior.pois2,
      #random =~ us(trait):mm(Sp + Sp2),
      family = "poisson",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf), file = paste0("Model Runs/Poisson Model ",i, ".Rdata"))
  }
}, mc.cores = 20)