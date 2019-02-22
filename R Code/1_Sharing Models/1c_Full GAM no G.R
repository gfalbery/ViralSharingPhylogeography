
# Fitting STAN Model in BRMS
# Rscript "R Code/1_Sharing Models/1c_Full GAM No G.R"

#source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew); library(brms)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables sans random effects

# Code for BRMS ####

BinGAMNoG <- brm(VirusBinary ~ t2(Space, scale(Phylo2)) + s(DietSim), # + 
                #(1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
              data = FinalHostMatrix, 
              family = bernoulli(), 
              cores = 4,
              iter = 750, 
              warmup = 250, 
              thin = 10, 
              refresh = 0,
              silent = FALSE, 
              verbose = TRUE)

saveRDS(BinGAM, file = "BinGAMNoG.rds")

library(mgcv)
