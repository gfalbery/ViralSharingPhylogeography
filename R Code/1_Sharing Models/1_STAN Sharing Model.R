
# STAN Model ####

# nice -n 10 Rscript "R Code/1_Sharing Models/STAN Model.R" # This is the terminal run code

source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

d <- FinalHostMatrix

# Do we have complete data for all rows?

sum(complete.cases(d)) == nrow(d)

# Generate species-level trait data

species.traits <- 
  data.frame(sp = c(d$Sp, d$Sp2),
             d_cites = c(log(d$hDiseaseZACites+1),log(d$hDiseaseZACites.Sp2+1)),
             domestic = c(d$hDom, d$hDom.Sp2)) %>%
  arrange(sp) %>%
  distinct() %>%
  mutate(sp = as.factor(sp),
         domestic = ifelse(domestic == "domestic", 1, 0))

# Get Sp and Sp2 in "d" on the same factor levels

d$Sp <- factor(as.character(d$Sp),
               levels = union(d$Sp, d$Sp2)
)

d$Sp2 <- factor(as.character(d$Sp2),
                levels = union(d$Sp, d$Sp2)
)

# Generate Stan data
stan.data <- list(
  
  N = nrow(d),
  
  virus_shared = d$VirusBinary,
  
  space_s = c(scale(d$Space)),
  
  phylo_s = c(scale(d$Phylo2)),
  
  #domdom = d$DomDom,
  
  d_cites_s = c(scale(species.traits$d_cites)),
  
  domestic = species.traits$domestic,
  
  species1 = to_stan_factor(d$Sp),
  species2 = to_stan_factor(d$Sp2),
  
  N_species = stan_factor_count(d$Sp2)
) %>% plyr::mutate(space_phylo_s = c(scale(d$Space*d$Phylo2)))

# Set Stan model parameters

iter <- 1500
warmup <- 500
chains <- 8
cores <- 8
adapt_delta <- 0.99
stepsize <- 0.5

BinModel_Scaled <- 
  stan(file = "R Code/1_Sharing Models/Albersnet.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(BinModel_Scaled, 
        file = "BinModel_Scaled.rds")

