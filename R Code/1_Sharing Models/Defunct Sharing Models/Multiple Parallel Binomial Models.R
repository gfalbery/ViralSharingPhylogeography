# Parallel binary models ####

#nice -n 10 Rscript "R Code/1_Sharing Models/Multiple Parallel Binomial Models.R" # This is the terminal run code

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

BinModelList <- BinModelList2 <- BinModelList3 <- list()

BinModelList <- parallel::mclapply(1:20, function(i) {
  if(i <= 10) {
#  
#    MCMCglmm(
#      data = FinalHostMatrix,
#     VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
#      prior = prior.bin,
#      random =~ mm(Sp + Sp2),
#      family = "categorical",
#      pr = TRUE,
#      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
#      thin = 10*mf, burnin=8000*mf, trunc = T) %>% return
    
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

save(BinModelList, file = "Parallel_Binomialsb.Rdata")

BinModelList2 <- parallel::mclapply(21:40, function(i) {
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

save(BinModelList2, file = "Parallel_Binomials2.Rdata")

BinModelList3 <- parallel::mclapply(41:50, function(i) {
  
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

save(BinModelList3, file = "Parallel_Binomials3.Rdata")
