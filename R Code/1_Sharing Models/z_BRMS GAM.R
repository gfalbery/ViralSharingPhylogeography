
# Fitting STAN Model in BRMS
# Rscript "R Code/1_Sharing Models/1b_Full GAM.R"
# Rscript --verbose "R Code/1_Sharing Models/z_BRMS GAM.R" >> my_log.txt

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew); library(brms)

# Binomial model for viral sharing with varying effects by species 
# modeled using predictor variables

# Import data

f <- FinalHostMatrix

# Get Sp and Sp2 in "d" on the same factor levels

f$Sp <- factor(as.character(f$Sp),
               levels = union(f$Sp, f$Sp2)
)

f$Sp2 <- factor(as.character(f$Sp2),
                levels = union(f$Sp, f$Sp2)
)

# Generate Stan data
f <- f %>% mutate(
  
  space_s = f$Space,
  phylo_s = c(scale(f$Phylo2)),
  
  dom = ifelse(hDom=="wild",0,1),
  dom_2 = ifelse(hDom.Sp2=="wild",0,1),
  
  d_cites_s1 = scale(log(hDiseaseZACites+1))[,1],
  d_cites_s2 = scale(log(hDiseaseZACites.Sp2+1))[,1],
  
)

# Code for BRMS ####

sink(stdout())

BRMSGAM <- brm(VirusBinary ~ t2(space_s,phylo_s) + s(DietSim) +
              (1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
            data = f, 
            family = bernoulli(), 
            cores = 10,
            chains = 10,
            iter = 1000,
            warmup = 500,
            thin = 10, 
            refresh = 0,
            save_ranef = TRUE,
            verbose = TRUE,
            silent = FALSE)

sink()

save(BRMSGAM, file = "Model Files/BRMSGAM.Rdata")


