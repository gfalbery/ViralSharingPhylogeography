# Extracting deviance parameters from BAM ####

# Rscript "R Code/1_Sharing Models/1c_Frequentist GAM Deviance Constributions.R" BAM ####

# Running Frequentist GAMS

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

load("Output Files/BAMList.Rdata")
#load("Output Files/BAMList2.Rdata")

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

library(mgcv); library(tidyverse)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

DevList <- DevList2 <- DataList <- PPList <- list()

for(r in 1:length(BAMList)){
  
  print(Resps[r])
  
  DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r]))
  
  DataList[[Resps[r]]]$Sp <- factor(DataList[[Resps[r]]]$Sp, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  DataList[[Resps[r]]]$Sp2 <- factor(DataList[[Resps[r]]]$Sp2, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  
  MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  
  SppMatrix = MZ1 + MZ2
  
  DataList[[Resps[[r]]]]$Spp <- SppMatrix
  DataList[[Resps[[r]]]]$Cites <- rowSums(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
  DataList[[Resps[[r]]]]$MinCites <- apply(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
  DataList[[Resps[[r]]]]$Domestic <- ifelse(rowSums(cbind(2- FinalHostMatrix$hDom %>% as.factor %>% as.numeric,
                                                          2- FinalHostMatrix$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)
  
  PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), diag(nlevels(DataList[[Resps[r]]]$Sp))))
}

Covar <- c("s(Space)", 
           "s(scale(Phylo2))",
           "s(DietSim)",
           "MinCites", 
           "Domestic",
           "Spp")

Formula = as.formula(paste0(Resps[1], 
                            " ~ ",
                            paste(Covar, collapse = " + ")))

DevList$FullModel <- bam(Formula,
                         data = DataList[[Resps[1]]], 
                         family = binomial(),
                         paraPen = PPList[[Resps[1]]])

r = 1

for(r in 1:length(Covar)){
  
  print(Covar[r])
  
  TestCovar <- setdiff(Covar, Covar[r])
  
  Formula = as.formula(paste0(Resps[1], 
                              " ~ ",
                              paste(TestCovar, collapse = " + ")))
  
  DevList[[Covar[r]]] <- bam(Formula,
                             data = DataList[[Resps[1]]], 
                             family = binomial(),
                             paraPen = PPList[[Resps[1]]])
  
}

save(DevList, file = "Output Files/DevList.Rdata")

OrigDev = deviance(DevList$FullModel)
RemoveDev = sapply(DevList[2:length(DevList)], deviance)

DevExplained = (RemoveDev - OrigDev)

DevExplained/sum(DevExplained)



