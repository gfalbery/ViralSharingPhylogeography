
# Running Frequentist GAMS

# Rscript "R Code/1_Sharing Models/1a_Frequentist GAMs.R"

#if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else{
  source("R Code/00_Master Code.R")
#}

library(mgcv); library(tidyverse); library(ggregplot)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

Covar <- c("s(Phylo, by = ordered(Gz))",
           "t2(Phylo, Space, by = ordered(!Gz))",
           "MinCites",
           "Domestic",
           "Spp")

BAMList <- DataList <- PPList <- list()

r = 1

for(r in 1:length(Resps)){
  
  print(Resps[r])
  
  DataList[[Resps[r]]] <- FinalHostMatrix[!NARows(FinalHostMatrix, Resps[r]),] %>% droplevels
  
  DataList[[Resps[r]]]$Sp <- factor(DataList[[Resps[r]]]$Sp, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  DataList[[Resps[r]]]$Sp2 <- factor(DataList[[Resps[r]]]$Sp2, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  
  DataList[[Resps[r]]] <- DataList[[Resps[r]]] %>% slice(order(Sp, Sp2))
  
  MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  
  SppMatrix = MZ1 + MZ2
  
  DataList[[Resps[[r]]]]$Spp <- SppMatrix
  DataList[[Resps[[r]]]]$Cites <- rowSums(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
  DataList[[Resps[[r]]]]$MinCites <- apply(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
  DataList[[Resps[[r]]]]$Domestic <- ifelse(rowSums(cbind(2- DataList[[Resps[r]]]$hDom %>% as.factor %>% as.numeric,
                                                          2- DataList[[Resps[r]]]$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)
  
  DataList[[Resps[[r]]]]$Space <- DataList[[Resps[[r]]]]$SpaceMin
  # DataList[[Resps[[r]]]]$Space <- DataList[[Resps[[r]]]]$SpaceMax
  
  PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), 
                                        diag(nlevels(DataList[[Resps[r]]]$Sp))))

  Formula = as.formula(paste0(Resps[r], 
                              " ~ ",
                              paste(Covar, collapse = " + ")
  ))
  
  BAMList[[Resps[r]]] <- bam(Formula,
                             data = DataList[[Resps[r]]], 
                             family = binomial(),
                             paraPen = PPList[[Resps[r]]], select = T
  )
  
}

save(DataList, PPList, BAMList, file = "Output Files/BAMList.Rdata")
