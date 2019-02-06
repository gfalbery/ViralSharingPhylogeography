# Parallel binary models ####

source("R Code/00_Master Code.R")

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel); library(dplyr)

prior.bin <- list(R = list(V = diag(1), nu = 0.002, fix = 1),
                  G = list(G1 = list(V = diag(1), nu = 2)))

prior.bin2 <- list(R = list(V = diag(1), nu = 0.002, fix = 1))

# Modelling all mammal-mammal pairs ####

mf = 1

# Trying a Binomial model ####

mc1 <- MCMCglmm(
  data = FinalHostMatrix,
  VirusBinary ~ 1,#Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
  random =~ mm(hDiseaseZACites + hDiseaseZACites.Sp2)*mm(Sp + Sp2),
  prior = prior.bin,
  family = "categorical",
  pr = TRUE,
  nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin=3000*mf, trunc = T)

save(mc1, file = "Bin Model 1.Rdata")

mc2 <- MCMCglmm(
  data = FinalHostMatrix,
  VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
  prior = prior.bin2,
  family = "categorical",
  #pr = TRUE,
  nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin=3000*mf, trunc = T)

mc2b <- MCMCglmm(
  data = FinalHostMatrix,
  VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
  prior = prior.bin2,
  family = "categorical",
  #pr = TRUE,
  nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin=3000*mf, trunc = F)

save(mc2, file = "Bin Model 2.Rdata")


