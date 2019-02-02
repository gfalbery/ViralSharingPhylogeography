# Parallel binary models ####

# Run source code

rm(list = ls())

source("R Code/00_Master Code.R")

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel); library(tidyverse)

prior.bin <- list(R = list(V = diag(1), nu = 0.002, fix = 1),
                  G = list(G1 = list(V = diag(1), nu = 2)))

prior.bin2 <- list(R = list(V = diag(1), nu = 0.002, fix = 1))

# Modelling all mammal-mammal pairs ####

mf = 15

# Trying a Binomial model ####

BinModelList <- list()

BinModelList[1:10] <- NA

BinModelList <- parallel::mclapply(1:20, function(i) {
  if(i <= 10) {
  
    MCMCglmm(
      data = FinalHostMatrix,
     VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      prior = prior.bin,
      random =~ mm(Sp + Sp2),
      family = "categorical",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, trunc = T) %>% return
    
  } else if (i > 10) {
    
    MCMCglmm(
      data = FinalHostMatrix,
      VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      prior = prior.bin2,
      family = "categorical",
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, trunc = T) %>% return
  }
}, mc.cores = 20)

save(BinModelList, file = "Parallel_Binomials.Rdata")

BinModelList[21:40] <- parallel::mclapply(21:40, function(i) {
  if(i <= 30) {
    
    MCMCglmm(
      data = FinalHostMatrixNoSpace,
      VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      prior = prior.bin,
      random =~ mm(Sp + Sp2),
      family = "categorical",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, trunc = T) %>% return
    
  } else if (i > 30) {
    
    MCMCglmm(
      data = FinalHostMatrixNoSpace,
      VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      prior = prior.bin2,
      family = "categorical",
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, trunc = T) %>% return
    
  }
}, mc.cores = 20)

save(BinModelList[21:40], file = "Parallel_Binomials2.Rdata")

BinModelList[41:50] <- parallel::mclapply(41:50, function(i) {
  
  MCMCglmm(
    data = FinalHostMatrix,
    VirusBinary ~ MinDCites,
    prior = prior.bin,
    random =~ mm(Sp + Sp2),
    family = "categorical",
    pr = TRUE,
    nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
    thin = 10*mf, burnin=8000*mf, trunc = T) %>% return
  
}, mc.cores = 10)

save(BinModelList[41:50], file = "Parallel_Binomials3.Rdata")
