
# Fitting STAN Model in BRMS

source("R Code/00_Master Code.R")

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

m3.3 <- brm(DNA ~ t2(space_s,phylo_s) + 
              (1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
            data = f, 
            family = bernoulli(), 
            cores = 8, 
            seed = 17,
            iter = 1000, 
            warmup = 250, 
            thin = 10, 
            refresh = 0)

saveRDS(m3.3, file = "DNAGAM.rds")

VectorGAM <- brm(Vector ~ t2(space_s,phylo_s) + 
                   (1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
                 data = FinalHostMatrix %>% filter(!is.na(Vector)), 
                 family = bernoulli(), 
                 cores = 8, 
                 seed = 17,
                 iter = 1000, 
                 warmup = 250, 
                 thin = 10, 
                 refresh = 0)

saveRDS(VectorGAM, file = "VectorGAM.rds")

NVectorGAM <- brm(NVector ~ t2(space_s,phylo_s) + 
                    (1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
                  data = FinalHostMatrix %>% filter(!is.na(NVector)), 
                  family = bernoulli(), 
                  cores = 8, 
                  seed = 17,
                  iter = 1000, 
                  warmup = 250, 
                  thin = 10, 
                  refresh = 0)

saveRDS(NVectorGAM, file = "NVectorGAM.rds")

plot(marginal_smooths(m3.3), ask = FALSE)

