
# STAN Model ####

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
             d_cites = log(c(d$hDiseaseZACites, d$hDiseaseZACites.Sp2)+1),
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
assert_that(stan.data$N_species == n_species)

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

binom.ranef2.fit.model <- 
  stan(file = "albers.binom.ranef2.stan",
       data = stan.data,
       iter = iter, warmup = warmup,
       chains = chains, cores = cores, 
       control = list(
         adapt_delta = adapt_delta, stepsize = stepsize
       )
  )

saveRDS(binom.ranef2.fit.model, 
        file = "stan/saved_models/binom.ranef2.fit.model.rds")


traceplot(binom.ranef2.fit.model)
traceplot(binom.ranef2.fit.model,
          pars = c("mu_alpha", 
                   "beta_d_cites_s", 
                   "beta_domestic", 
                   "beta_space",
                   "beta_phylo",
                   "beta_inter",
                   "sigma"))

p <- process_stanfit(binom.ranef2.fit.model, n.pars.to.trim = 3)

precis(p$df, prob = 0.95)

species.means <- apply(p$df, 2, mean)
species.means <- species.means[c(-1, -651, -652, -653)]

plot(species.traits$d_cites_standardized, species.means,
     pch = 19, col = alpha("grey", 0.8))
plot(species.traits$domestic, species.means,
     pch = 19, col = alpha("grey", 0.8))


# Predictions!

# Packge together species-level varying effect estimates

species.means.df <- data.frame(
  species_mean = species.means,
  species_level = 1:length(species.means)
)

species.traits.mod <- species.traits %>%
  mutate(sp = as.integer(sp))

# Create a predictions data frame

preds <- data.frame(
  sp = d$Sp,
  sp_level = as.integer(d$Sp),
  sp2 = d$Sp2,
  sp2_level = as.integer(d$Sp2),
  # Mean overall intercept fit from the model
  mu_alpha = rep(mean(p$df$mu_alpha), nrow(d))
) %>%
  left_join(., species.means.df, by = c("sp_level" = "species_level")) %>%
  left_join(., species.means.df, by = c("sp2_level" = "species_level")) %>%
  left_join(., species.traits.mod, by = c("sp_level" = "sp")) %>%
  left_join(., species.traits.mod, by = c("sp2_level" = "sp")) %>%
  rename(
    # Mean species-level varying effects fit from the model
    fit_alpha_sp1 = species_mean.x,
    fit_alpha_sp2 = species_mean.y,
    # Species trait data
    domestic_sp1 = domestic.x,
    domestic_sp2 = domestic.y,
    d_cites_standardized_sp1 = d_cites_standardized.x,
    d_cites_standardized_sp2 = d_cites_standardized.y
  ) %>%
  mutate(
    # Manually calculate linear term that feeds into the varying effect
    # for each species. Note: will not match the fitted values, since the
    # fitted values are also subject to residual variation
    linear_term_alpha_sp1 =
      (domestic_sp1 * mean(p$df$beta_domestic) + 
         d_cites_standardized_sp1 * mean(p$df$beta_d_cites_s)),
    linear_term_alpha_sp2 =
      (domestic_sp2 * mean(p$df$beta_domestic) + 
         d_cites_standardized_sp2 * mean(p$df$beta_d_cites_s)),
    # Generate predicted varying effect values for species just like ours
    # (including the residual variation)
    pred_alpha_sp1 = rnorm(linear_term_alpha_sp1, sd = p$df$sigma),
    pred_alpha_sp2 = rnorm(linear_term_alpha_sp2, sd = p$df$sigma)
  )

# Generate predictions

# Predictions using mean fitted varying effect values from the model
preds$pred <- 
  rbinom(nrow(preds), 
         prob = logistic(preds$mu_alpha + preds$fit_alpha_sp1 + preds$fit_alpha_sp2), 
         size = 1
  )

# Predictions using manually calculated linear terms for our species
# (not subject to residual variation)
preds$pred2 <-
  rbinom(nrow(preds),
         prob = logistic(preds$mu_alpha + preds$linear_term_alpha_sp1 + preds$linear_term_alpha_sp2),
         size = 1
  )

# Predictions using manually calculated linear terms for our species
# (and the influence of residual variation)
preds$pred3 <-
  rbinom(nrow(preds),
         prob = logistic(preds$mu_alpha + preds$pred_alpha_sp1 + preds$pred_alpha_sp2),
         size = 1
  )


# How many viral sharing links are observed?

sum(d$VirusBinary)/nrow(d)

sum(preds$pred)/nrow(preds)
sum(preds$pred2)/nrow(preds)
sum(preds$pred3)/nrow(preds)