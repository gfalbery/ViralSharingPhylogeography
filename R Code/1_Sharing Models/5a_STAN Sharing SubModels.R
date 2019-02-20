
# STAN Model ####

# nice -n 10 Rscript "R Code/1_Sharing Models/5a_Sharing SubModels.R" # This is the terminal run code

source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew)

SubResps <- c("RNA","Vector","NVector","DNA")

# Binomial model for viral sharing of viral subtypes ####

SubDataList <- StanDataList <- SubModelList <- list()

# Import data

for(r in 1:length(SubResps)){
  
  SubDataList[[r]] <- FinalHostMatrix %>% filter(!is.na(SubResps[r]))
  
  # Generate species-level trait data
  
  species.traits <- 
    data.frame(sp = c(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2),
               d_cites = c(log(SubDataList[[r]]$hDiseaseZACites+1),log(SubDataList[[r]]$hDiseaseZACites.Sp2+1)),
               domestic = c(SubDataList[[r]]$hDom, SubDataList[[r]]$hDom.Sp2)) %>%
    arrange(sp) %>%
    distinct() %>%
    mutate(sp = as.factor(sp),
           domestic = ifelse(domestic == "domestic", 1, 0))
  
  # Get Sp and Sp2 in "d" on the same factor levels
  
  SubDataList[[r]]$Sp <- factor(as.character(SubDataList[[r]]$Sp),
                                levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  )
  
  SubDataList[[r]]$Sp2 <- factor(as.character(SubDataList[[r]]$Sp2),
                                 levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  )
  
  # Generate Stan data
  StanDataList[[r]] <- list(
    
    N = nrow(SubDataList[[r]]),
    
    virus_shared = SubDataList[[r]][,Resps[r]],
    
    space_s = c(scale(SubDataList[[r]]$Space)),
    
    phylo_s = c(scale(SubDataList[[r]]$Phylo2)),
    
    #domdom = SubDataList[[r]]$DomDom,
    
    d_cites_s = c(scale(species.traits$d_cites)),
    
    domestic = species.traits$domestic,
    
    species1 = to_stan_factor(SubDataList[[r]]$Sp),
    species2 = to_stan_factor(SubDataList[[r]]$Sp2),
    
    N_species = stan_factor_count(SubDataList[[r]]$Sp2)
  ) %>% plyr::mutate(space_phylo_s = c(scale(SubDataList[[r]]$Space*SubDataList[[r]]$Phylo2)))
  
  # Set Stan model parameters
  
  iter <- 1500
  warmup <- 500
  chains <- 8
  cores <- 8
  adapt_delta <- 0.99
  stepsize <- 0.5
  
  SubModel <- 
    stan(file = "R Code/1_Sharing Models/Albersnet.stan",
         data = stan.data,
         iter = iter, warmup = warmup,
         chains = chains, cores = cores, 
         control = list(
           adapt_delta = adapt_delta, stepsize = stepsize
         )
    )
  
  SubModelList[[r]] <- SubModel
  
  saveRDS(SubModel, 
          file = paste0(Resps[r],"SubModel.rds"))
  
  remove(SubModel)
  
}

saveRDS(SubModelList, file = "SubModelList.rds")
saveRDS(StanDataList, file = "StanDataList.rds")
saveRDS(SubDataList, file = "SubDataList.rds")