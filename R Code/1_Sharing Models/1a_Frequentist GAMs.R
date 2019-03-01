
# Running Frequentist GAMS

# Rscript "R Code/1_Sharing Models/1a_Frequentist GAMs.R"

#if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else 
source("R Code/00_Master Code.R")

library(mgcv); library(tidyverse)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

BAMList <- BAMList2 <- DataList <- PPList <- list()

for(r in 1:length(Resps)){
  
  print(Resps[r])
  
  DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r]))
  
  MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  
  SppMatrix = MZ1 + MZ2
  
  DataList[[Resps[[r]]]]$Spp <- SppMatrix
  DataList[[Resps[[r]]]]$Cites <- rowSums(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
  DataList[[Resps[[r]]]]$MinCites <- apply(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
  DataList[[Resps[[r]]]]$Domestic <- ifelse(rowSums(cbind(2- FinalHostMatrix$hDom %>% as.factor %>% as.numeric,
                                                          2- FinalHostMatrix$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)
  
  PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), 
                                        diag(nlevels(DataList[[Resps[r]]]$Sp))))
  
  Formula = as.formula(paste0(Resps[r], 
                              " ~ ",
                              #"t2(Space, scale(Phylo2)) + 
                              "s(Space) + s(scale(Phylo2)) + s(DietSim) + 
                              MinCites + Domestic + 
                              Spp"))
  
  Formula2 = as.formula(paste0(Resps[r], 
                               " ~ ",
                               #t2(Space, scale(Phylo2)) + 
                               "s(Space) + s(scale(Phylo2)) + s(DietSim) + 
                               MinCites + Domestic"))
  
  BAMList[[Resps[r]]] <- bam(Formula,
                             data = DataList[[Resps[r]]], 
                             family = binomial(),
                             paraPen = PPList[[Resps[r]]])
  
  BAMList2[[Resps[r]]] <- bam(Formula2,
                              data = DataList[[Resps[r]]], 
                              family = binomial())
  
}

save(BAMList, file = "Output Files/BAMList.Rdata")
save(BAMList2, file = "Output Files/BAMList2.Rdata")

