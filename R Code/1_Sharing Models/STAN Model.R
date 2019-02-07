
# STAN Model ####

# nice -n 10 Rscript "STAN Model.R" # This is the terminal run code

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
         #d_cites_standardized = scale(d_cites)[ ,1],
         domestic = ifelse(domestic == "domestic", 1, 0))

# Get Sp and Sp2 in "d" on the same factor levels

d$Sp <- factor(as.character(d$Sp),
               levels = union(d$Sp, d$Sp2)
)

d$Sp2 <- factor(as.character(d$Sp2),
                levels = union(d$Sp, d$Sp2)
)

summary(as.integer(d$Sp))
summary(as.integer(d$Sp2))

# Do the factors in Sp and Sp2 match with what's in species.traits?

n_species <- length(levels(d$Sp))

sum(levels(d$Sp) == levels(species.traits$sp)) == n_species
sum(levels(d$Sp2) == levels(species.traits$sp)) == n_species

# Generate Stan data
stan.data <- list(
  
  N = nrow(d),
  
  virus_shared = d$VirusBinary,
  
  space = d$Space,
  
  phylo = d$Phylo2,
  
  space_phylo = d$Space*d$Phylo2,
  
  #domdom = d$DomDom,
  
  d_cites_s = species.traits$d_cites,
  
  domestic = species.traits$domestic,
  
  species1 = to_stan_factor(d$Sp),
  
  species2 = to_stan_factor(d$Sp2),
  
  N_species = stan_factor_count(d$Sp2)
)

# Check to make sure the Stan factor variables are consistent across the
# two species columns
d$Sp2[1]
d$Sp[649] # Same as d$Sp2[1]

stan.data$species2[1]
stan.data$species1[649] # Should be same as stan.data$species2[1]

# Set Stan model parameters

iter <- 2500
warmup <- 1500
chains <- 8
cores <- 8
adapt_delta <- 0.99
stepsize <- 0.5

binom.model <- 
  stan(file = "Albersnet.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(binom.model, 
        file = "Bin Model.rds")

