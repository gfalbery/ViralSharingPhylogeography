
# Fitting STAN Model in BRMS
# Rscript "R Code/1_Sharing Models/1b_Full GAM.R"
# Rscript --verbose "R Code/1_Sharing Models/z_BRMS GAM.R" >> my_log.txt

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew); library(brms)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

load("Output Files/BAMList.Rdata")

BRMSData <- DataList[[1]] %>%
  mutate(DCites1 = log(hDiseaseZACites + 1),
         DCites2 = log(hDiseaseZACites.Sp2 + 1),
         Dom1 = ifelse(hDom == "domestic", 1, 0),
         Dom2 = ifelse(hDom.Sp2 == "domestic", 1, 0))

# Code for BRMS ####

sink(stdout())

BRMSGAM <- brm(VirusBinary ~ t2(Space, Phylo) + s(DietSim) +
                 (1 + mmc(Dom1, Dom2) + mmc(DCites1, DCites2) | mm(Sp, Sp2)),
               data = BRMSData, 
               family = bernoulli(), 
               cores = 10,
               chains = 10,
               iter = 2000,
               warmup = 750,
               thin = 10, 
               refresh = 0,
               save_ranef = TRUE,
               verbose = TRUE,
               silent = FALSE)

sink()

save(BRMSGAM, file = "Model Files/BRMSGAM.Rdata")


