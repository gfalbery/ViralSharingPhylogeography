
# Running RNA-specific model ####

# STAN Model ####

# Rscript "R Code/1_Sharing Models/1e_DNA STAN Model.R" # This is the terminal run code

source("R Code/00_Master Code.R")
source("R Code/1_Sharing Models/1d_Separating RNA and DNA.R")

library(rstan); library(tidyverse); library(reskew)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

f <- FinalHostMatrix %>% filter(!is.na(DNA))

# Generate species-level trait data

species.traits <- 
  data.frame(sp = c(f$Sp, f$Sp2),
             d_cites = c(log(f$hDiseaseZACites+1),log(f$hDiseaseZACites.Sp2+1)),
             domestic = c(f$hDom, f$hDom.Sp2)) %>%
  arrange(sp) %>%
  distinct() %>%
  mutate(sp = as.factor(sp),
         domestic = ifelse(domestic == "domestic", 1, 0))

# Get Sp and Sp2 in "d" on the same factor levels

f$Sp <- factor(as.character(f$Sp),
               levels = union(f$Sp, f$Sp2)
)

f$Sp2 <- factor(as.character(f$Sp2),
                levels = union(f$Sp, f$Sp2)
)

# Generate Stan data
stan.data <- list(
  
  N = nrow(f),
  
  DNAvirus_shared = f$DNA,
  
  space_s = c(scale(f$Space)),
  
  phylo_s = c(scale(f$Phylo2)),
  
  #domdom = d$DomDom,
  
  d_cites_s = c(scale(species.traits$d_cites)),
  
  domestic = species.traits$domestic,
  
  species1 = to_stan_factor(f$Sp),
  species2 = to_stan_factor(f$Sp2),
  
  N_species = stan_factor_count(f$Sp2)
) %>% plyr::mutate(space_phylo_s = c(scale(f$Space*f$Phylo2)))

# Set Stan model parameters

iter <- 1000
warmup <- 250
chains <- 8
cores <- 8
adapt_delta <- 0.99
stepsize <- 0.5

DNABinModel <- 
  stan(file = "R Code/1_Sharing Models/DNABinModel.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(DNABinModel, 
        file = "DNABinModel.rds")

