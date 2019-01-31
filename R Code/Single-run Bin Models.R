# Parallel binary models ####

# Run parallel

library(MCMCglmm); library(ggregplot); library(INLA); library(parallel); library(dplyr)

prior.bin2 <- list(R = list(V = diag(1), nu = 0.002, fix = 1))

# Modelling all mammal-mammal pairs ####

mf = 10

# Trying a Binomial model ####

mc2 <- MCMCglmm(
  data = FinalHostMatrix,
  VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
  prior = prior.bin2,
  family = "categorical",
  pr = TRUE,
  nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin=3000*mf, trunc = T)


mc1 <- MCMCglmm(
  data = FinalHostMatrix,
  VirusBinary ~ Space + Phylo2 + Space:Phylo2 + MinDCites + DomDom,
  random =~ mm(Sp + Sp2),
  prior = prior.bin,
  family = "categorical",
  pr = TRUE,
  nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin=3000*mf, trunc = T)

