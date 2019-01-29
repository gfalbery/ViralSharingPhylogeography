
# Shorter parallel
# Because the spatial overlap model won't run without Sp and Sp2 having the same levels
# Run source code

rm(list = ls())

source("R Code/00_Master Code.R")

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel); library(dplyr)

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo2
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrixNoSpace <- FinalHostMatrix %>% filter(Space>0)

OverlapMatrix <- droplevels(FinalHostMatrix[FinalHostMatrix$Space>0,])
OverlapMatrix$Sp <- factor(OverlapMatrix$Sp, levels = union(OverlapMatrix$Sp,OverlapMatrix$Sp2))
OverlapMatrix$Sp2 <- factor(OverlapMatrix$Sp2, levels = union(OverlapMatrix$Sp,OverlapMatrix$Sp2))

gp <- gelman.prior(Virus ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom, data = FinalHostMatrix)

prior.zi <- list(R = list(V = diag(2), nu = 0, fix = 2),
                 G = list(G1 = list(V = diag(2), nu = 2)),
                 B = list(V = diag(14)*10^8, mu = rep(0,14)))

prior.zi$B$V[seq(2,dim(gp)[2]*2,2),seq(2,dim(gp)[2]*2,2)] <- gp

mf = 15

ZI_shortruns <- parallel::mclapply(1:10, function(i) {
  
  return(MCMCglmm(
    data = OverlapMatrix,
    Virus ~ trait -1 + trait:(Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom),
    rcov =~ idh(trait):units, 
    prior = prior.zi,
    random =~ us(trait):mm(Sp + Sp2),
    family = "zipoisson",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf))}, mc.cores = 10)

# Save file

save(ZI_shortruns, file = "Model Runs/ZI_shortruns.Rdata")
