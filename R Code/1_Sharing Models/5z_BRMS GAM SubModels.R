
# GAM SubModels ####

# Rscript "~/Albersnet/R Code/1_Sharing Models/5z_BRMS GAM SubModels.R" # This is the terminal run code

source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew); library(brms)

SubResps <- c("RNA","Vector","NVector","DNA")

# Binomial model for viral sharing of viral subtypes ####

SubDataList <- StanDataList <- SubGAMModelList <- list()

# Import data

for(r in 1:length(SubResps[1:3])){
  
  SubDataList[[r]] <- FinalHostMatrix[!is.na(FinalHostMatrix[,SubResps[r]]),]
  
  # Generate species-level trait data
    # Get Sp and Sp2 in "d" on the same factor levels
  
  SubDataList[[r]]$Sp <- factor(as.character(SubDataList[[r]]$Sp),
                                levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  )
  
  SubDataList[[r]]$Sp2 <- factor(as.character(SubDataList[[r]]$Sp2),
                                 levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  )
  
  SubDataList[[r]]$Sharing <- SubDataList[[r]][,SubResps[r]]
  
  SubDataList[[r]]$space_s <- scale(SubDataList[[r]]$Space) %>% c
  SubDataList[[r]]$phylo_s <- scale(SubDataList[[r]]$Phylo2) %>% c
  SubDataList[[r]]$d_cites_s <- scale(log(SubDataList[[r]]$hDiseaseZACites+1)) %>% c
  SubDataList[[r]]$d_cites_s2 <- scale(log(SubDataList[[r]]$hDiseaseZACites.Sp2+1)) %>% c
  
  SubDataList[[r]]$domestic <- ifelse(SubDataList[[r]]$hDom=="domestic",1,0)
  SubDataList[[r]]$domestic.Sp2 <- ifelse(SubDataList[[r]]$hDom.Sp2=="domestic",1,0)
  
  SubModel <- brm(Sharing ~ t2(space_s, phylo_s) + 
                (1 + mmc(domestic, domestic.Sp2) + mmc(d_cites_s, d_cites_s2) | mm(Sp, Sp2)),
                data = SubDataList[[r]], 
                family = bernoulli(), 
                cores = 8, 
                seed = 17,
                iter = 1500, 
                warmup = 500, 
                thin = 10, 
                refresh = 0)
  
  SubGAMModelList[[r]] <- SubModel
  
  saveRDS(SubModel, 
          file = paste0(Resps[r],"SubModel.rds"))
  
  remove(SubModel)
  
}

saveRDS(SubGAMModelList, file = "Output Files/SubModelList.rds")
saveRDS(SubDataList, file = "Output Files/SubDataList.rds")

