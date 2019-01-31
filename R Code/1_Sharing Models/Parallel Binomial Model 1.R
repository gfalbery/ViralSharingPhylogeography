# Parallel binary models ####

# Run source code

rm(list = ls())

source("R Code/00_Master Code.R")

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel); library(dplyr)

prior.bin <- list(R = list(V = diag(1), nu = 0.002, fix = 1),
                  G = list(G1 = list(V = diag(1), nu = 2)))

# Modelling all mammal-mammal pairs ####

mf = 15

# Trying a Binomial model ####

parallel::mclapply(1:10, function(i) {

    saveRDS(MCMCglmm(
      data = FinalHostMatrix,
      VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
      prior = prior.bin,
      random =~ mm(Sp + Sp2),
      family = "categorical",
      pr = TRUE,
      nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
      thin = 10*mf, burnin=8000*mf, trunc = T), file = paste0("Binomial Model ",i, ".Rdata"))
    
  }, mc.cores = 20)