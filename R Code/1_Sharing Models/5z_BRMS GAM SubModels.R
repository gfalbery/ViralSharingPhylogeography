
# GAM SubModels ####

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
    # Get Sp and Sp2 in "d" on the same factor levels
  
  SubDataList[[r]]$Sp <- factor(as.character(SubDataList[[r]]$Sp),
                                levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  )
  
  SubDataList[[r]]$Sp2 <- factor(as.character(SubDataList[[r]]$Sp2),
                                 levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  )
  
  SubDataList[[r]]$Sharing <- SubDataList[[r]][,SubResps[r]]
  
  SubModel <- brm(Sharing ~ t2(space_s, phylo_s), #+ 
                #(1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
                data = SubDataList[[r]], 
                family = bernoulli(), 
                cores = 8, 
                seed = 17,
                iter = 1000, 
                warmup = 250, 
                thin = 10, 
                refresh = 0)
  
  SubModelList[[r]] <- SubModel
  
  saveRDS(SubModel, 
          file = paste0(Resps[r],"SubModel.rds"))
  
  remove(SubModel)
  
}

saveRDS(SubModelList, file = "SubModelList.rds")
saveRDS(StanDataList, file = "StanDataList.rds")
saveRDS(SubDataList, file = "SubDataList.rds")

# Running it ####
