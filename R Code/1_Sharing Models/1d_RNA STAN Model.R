
# Running RNA-specific model ####

# STAN Model ####

# nice -n 10 Rscript "R Code/1_Sharing Models/1d_RNA STAN Model.R" # This is the terminal run code

source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

e <- FinalHostMatrix %>% filter(!is.na(RNA))

# Generate species-level trait data

species.traits <- 
  data.frame(sp = c(e$Sp, e$Sp2),
             d_cites = c(log(e$hDiseaseZACites+1),log(e$hDiseaseZACites.Sp2+1)),
             domestic = c(e$hDom, e$hDom.Sp2)) %>%
  arrange(sp) %>%
  distinct() %>%
  mutate(sp = as.factor(sp),
         domestic = ifelse(domestic == "domestic", 1, 0))

# Get Sp and Sp2 in "d" on the same factor levels

e$Sp <- factor(as.character(e$Sp),
               levels = union(e$Sp, e$Sp2)
)

e$Sp2 <- factor(as.character(e$Sp2),
                levels = union(e$Sp, e$Sp2)
)

# Generate Stan data
stan.data <- list(
  
  N = nrow(e),
  
  RNAvirus_shared = e$RNA,
  
  space_s = c(scale(e$Space)),
  
  phylo_s = c(scale(e$Phylo2)),
  
  #domdom = d$DomDom,
  
  d_cites_s = c(scale(species.traits$d_cites)),
  
  domestic = species.traits$domestic,
  
  species1 = to_stan_factor(e$Sp),
  species2 = to_stan_factor(e$Sp2),
  
  N_species = stan_factor_count(e$Sp2)
) %>% plyr::mutate(space_phylo_s = c(scale(e$Space*e$Phylo2)))

# Set Stan model parameters

iter <- 1500
warmup <- 500
chains <- 8
cores <- 8
adapt_delta <- 0.99
stepsize <- 0.5

RNABinModel <- 
  stan(file = "R Code/1_Sharing Models/RNABinModel.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(RNABinModel, 
        file = "RNABinModel.rds")

