# Extracting deviance parameters from BAM ####

# Rscript "R Code/1_Sharing Models/1c_Tensor Deviance Contributions.R"

# Running Frequentist GAMS

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

load("Output Files/BAMList.Rdata")
load("Output Files/BAMList2.Rdata")

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

library(mgcv); library(tidyverse)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")[1]

TensorDevList <- TensorDevList2 <- DataList <- PPList <- list()

for(r in 1:length(BAMList)){
  
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
  
  PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), diag(nlevels(DataList[[Resps[r]]]$Sp))))
}

Covar <- c("t2(Space, scale(Phylo2))", "s(Space)", "s(scale(Phylo2))","s(DietSim)",
           "MinCites", "Domestic","Spp")

r = 1

Formula <- as.formula(paste0(Resps[r], " ~ ", paste(Covar[c(1,4:7)], collapse = " + ")))

TensorDevList[["FullModel"]] <- bam(Formula,
                                    data = DataList[[Resps[r]]], 
                                    family = binomial(),
                                    paraPen = PPList[[Resps[r]]])

Formula <- as.formula(paste0(Resps[r], " ~ ", paste(Covar[2:7], collapse = " + ")))

TensorDevList[[Covar[1]]] <- bam(Formula,
                                 data = DataList[[Resps[r]]], 
                                 family = binomial(),
                                 paraPen = PPList[[Resps[r]]])

Covar2 <- Covar[c(2:7)]

for(r in 1:length(Covar2)){
  
  Covar3 <- Covar2
  
  if(r<3) Covar3 <- Covar2 else Covar3 <- Covar[c(1,4:7)]
  
  print(Covar2[r])
  
  TestCovar <- setdiff(Covar3, Covar2[r])
  
  Formula = as.formula(paste0(Resps[1], 
                              " ~ ",
                              paste(TestCovar, collapse = " + ")))
  
  TensorDevList[[Covar2[r]]] <- bam(Formula,
                                    data = DataList[[Resps[1]]], 
                                    family = binomial(),
                                    paraPen = PPList[[Resps[1]]])
  
}

save(TensorDevList, file = "Output Files/TensorDevList.Rdata")
