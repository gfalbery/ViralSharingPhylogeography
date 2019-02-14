
# Running RNA-specific model ####

# STAN Model ####

# Rscript "R Code/1_Sharing Models/1f_Vector STAN Model.R" # This is the terminal run code

source("R Code/00_Master Code.R")
source("R Code/0_Data Import/0j_Separating RNA and DNA.R")

library(rstan); library(tidyverse); library(reskew)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

g <- FinalHostMatrix %>% filter(!is.na(Vector))

# Generate species-level trait data

species.traits <- 
  data.frame(sp = c(g$Sp, g$Sp2),
             d_cites = c(log(g$hDiseaseZACites+1),log(g$hDiseaseZACites.Sp2+1)),
             domestic = c(g$hDom, g$hDom.Sp2)) %>%
  arrange(sp) %>%
  distinct() %>%
  mutate(sp = as.factor(sp),
         domestic = ifelse(domestic == "domestic", 1, 0))

# Get Sp and Sp2 in "d" on the same factor levels

g$Sp <- factor(as.character(g$Sp),
               levels = union(g$Sp, g$Sp2)
)

g$Sp2 <- factor(as.character(g$Sp2),
                levels = union(g$Sp, g$Sp2)
)

# Generate Stan data
stan.data <- list(
  
  N = nrow(g),
  
  DNAvirus_shared = g$Vector,
  
  space_s = c(scale(g$Space)),
  
  phylo_s = c(scale(g$Phylo2)),
  
  #domdom = d$DomDom,
  
  d_cites_s = c(scale(species.traits$d_cites)),
  
  domestic = species.traits$domestic,
  
  species1 = to_stan_factor(g$Sp),
  species2 = to_stan_factor(g$Sp2),
  
  N_species = stan_factor_count(g$Sp2)
) %>% plyr::mutate(space_phylo_s = c(scale(g$Space*g$Phylo2)))

# Set Stan model parameters

iter <- 1000
warmup <- 250
chains <- 8
cores <- 8
adapt_delta <- 0.99
stepsize <- 0.5

VectorBinModel <- 
  stan(file = "R Code/1_Sharing Models/DNABinModel.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(VectorBinModel, 
        file = "VectorBinModel.rds")

