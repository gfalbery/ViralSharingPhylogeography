
# Fitting STAN Model in BRMS


source("R Code/00_Master Code.R")
source("R Code/0_Data Import/0j_Separating RNA and DNA.R")

library(rstan); library(tidyverse); library(reskew); library(brms)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

f <- FinalHostMatrix %>% filter(!is.na(DNA))

# Get Sp and Sp2 in "d" on the same factor levels

f$Sp <- factor(as.character(f$Sp),
               levels = union(f$Sp, f$Sp2)
)

f$Sp2 <- factor(as.character(f$Sp2),
                levels = union(f$Sp, f$Sp2)
)

# Generate Stan data
f <- f %>% mutate(
  
  N = nrow(f),
  
  DNAvirus_shared = f$DNA,
  
  space_s = c(scale(f$Space)),
  phylo_s = c(scale(f$Phylo2)),
  
  dom = ifelse(hDom=="wild",0,1),
  dom_2 = ifelse(hDom.Sp2=="wild",0,1),
  
  d_cites_s1 = scale(hDiseaseZACites)[,1],
  d_cites_s2 = scale(hDiseaseZACites.Sp2)[,1],
  
) %>% 
  mutate(space_phylo_s = c(scale(f$Space*f$Phylo2)))


# Code for BRMS ####

m2 <- brm(DNA ~ space_s + phylo_s + space_phylo_s + 
            (1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
          data = f, 
          family = bernoulli(), 
          cores = 8, 
          seed = 17,
          iter = 1000, 
          warmup = 250, 
          thin = 10, 
          refresh = 0,
          adapt_delta = 0.95)


m3 <- lmer(DNA ~ space_s + phylo_s + space_phylo_s + 
             (1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
           data = f)

# for GAM
m2 <- brm(bf(accel ~ s(times)),
          data = mcycle, family = gaussian(), cores = 4, seed = 17,
          iter = 4000, warmup = 1000, thin = 10, refresh = 0,
          control = list(adapt_delta = 0.99))




