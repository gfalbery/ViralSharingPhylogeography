# Extracting deviance parameters from BAM ####

# Rscript "R Code/1_Sharing Models/1c_Tensor Deviance Contributions.R"

# Running Frequentist GAMS

#if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

load("Output Files/BAMList.Rdata")
#load("Output Files/BAMList2.Rdata")

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

library(mgcv); library(tidyverse)

TensorDevList <- DevOutput <- list()

for(r in 1:length(BAMList)){
  
  print(Resps[r])
  
  Covar <- c(#"ti(Space, Phylo)", 
    #"s(Space)", 
    #"s(Phylo)",
    "t2(Space, Phylo)",
    "s(DietSim)",
    "MinCites", "Domestic", "Spp")
  
  Formula <- as.formula(paste0(Resps[r], " ~ ", paste(Covar, collapse = " + ")))
  
  TensorDevList[[Resps[r]]][["FullModel"]] <- bam(Formula,
                                                  data = DataList[[Resps[r]]], 
                                                  family = binomial(),
                                                  paraPen = PPList[[Resps[r]]])
  
  for(s in 1:length(Covar)){
    
    print(Covar[s])
    
    TestCovar <- Covar[-s]
    
    Formula = as.formula(paste0(Resps[r], 
                                " ~ ",
                                paste(TestCovar, collapse = " + ")))
    
    TensorDevList[[Resps[r]]][[Covar[s]]] <- bam(Formula,
                                                 data = DataList[[Resps[r]]], 
                                                 family = binomial(),
                                                 paraPen = PPList[[Resps[r]]])
    
  }
  
  OrigDev <- TensorDevList[[Resps[r]]]$FullModel %>% deviance
  
  RemoveDevs <- sapply(TensorDevList[[Resps[r]]][2:length(TensorDevList[[Resps[r]]])], deviance)
  
  DevOutput[[Resps[r]]] <- 
    round(((RemoveDevs - OrigDev)/sum(RemoveDevs - OrigDev)), 3)
  
}

save(TensorDevList, DevOutput, file = "Output Files/TensorDevList.Rdata")

load("Output Files/TensorDevList.Rdata")

DevOutput <- list()
