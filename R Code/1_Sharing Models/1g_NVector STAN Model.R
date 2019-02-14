
# Running RNA-specific model ####

# STAN Model ####

# Rscript "R Code/1_Sharing Models/1g_NVector STAN Model.R" # This is the terminal run code

source("R Code/00_Master Code.R")
source("R Code/0_Data Import/0j_Separating RNA and DNA.R")

library(rstan); library(tidyverse); library(reskew)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

h <- FinalHostMatrix %>% filter(!is.na(NVector))

# Generate species-level trait data

species.traits <- 
  data.frame(sp = c(h$Sp, h$Sp2),
             d_cites = c(log(h$hDiseaseZACites+1),log(h$hDiseaseZACites.Sp2+1)),
             domestic = c(h$hDom, h$hDom.Sp2)) %>%
  arrange(sp) %>%
  distinct() %>%
  mutate(sp = as.factor(sp),
         domestic = ifelse(domestic == "domestic", 1, 0))

# Get Sp and Sp2 in "d" on the same factor levels

h$Sp <- factor(as.character(h$Sp),
               levels = union(h$Sp, h$Sp2)
)

h$Sp2 <- factor(as.character(h$Sp2),
                levels = union(h$Sp, h$Sp2)
)

# Generate Stan data
stan.data <- list(
  
  N = nrow(h),
  
  DNAvirus_shared = h$NVector,
  
  space_s = c(scale(h$Space)),
  
  phylo_s = c(scale(h$Phylo2)),
  
  #domdom = d$DomDom,
  
  d_cites_s = c(scale(species.traits$d_cites)),
  
  domestic = species.traits$domestic,
  
  species1 = to_stan_factor(h$Sp),
  species2 = to_stan_factor(h$Sp2),
  
  N_species = stan_factor_count(h$Sp2)
) %>% plyr::mutate(space_phylo_s = c(scale(h$Space*h$Phylo2)))

# Set Stan model parameters

iter <- 1000
warmup <- 250
chains <- 8
cores <- 8
adapt_delta <- 0.99
stepsize <- 0.5

NVectorBinModel <- 
  stan(file = "R Code/1_Sharing Models/DNABinModel.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(NVectorBinModel, 
        file = "NVectorBinModel.rds")

