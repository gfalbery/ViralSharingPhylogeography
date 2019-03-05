# Extracting deviance parameters from BAM ####

# Rscript "R Code/1_Sharing Models/1c_Tensor Deviance Contributions.R"

# Running Frequentist GAMS

#if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

load("Output Files/BAMList.Rdata")
load("Output Files/BAMList2.Rdata")

Resps = c("VirusBinary")

library(mgcv); library(tidyverse)

TensorDevList <- TensorDevList2 <- list()

for(r in 1){#:length(BAMList)){
  
  print(Resps[r])
  
  Covar <- c("t2(Space, Phylo)", "s(Space)", "s(Phylo)","s(DietSim)",
             "MinCites", "Domestic","Spp")
  
  Formula <- as.formula(paste0(Resps[r], " ~ ", paste(Covar[c(1,4:7)], collapse = " + ")))
  
  TensorDevList[["FullModel"]] <- bam(Formula,
                                      data = DataList[[Resps[r]]], 
                                      family = binomial(),
                                      paraPen = PPList[[Resps[r]]])
  
  Formula <- as.formula(paste0(Resps[r], " ~ ", paste(Covar[2:7], collapse = " + ")))
  
  TensorDevList[[Resps[r]]]$NoTensior <- bam(Formula,
                                             data = DataList[[Resps[r]]], 
                                             family = binomial(),
                                             paraPen = PPList[[Resps[r]]])
  
  Covar2 <- Covar[c(2:7)]
  
  for(s in 1:length(Covar2)){
    
    Covar3 <- Covar2
    
    if(s<3) Covar3 <- Covar2 else Covar3 <- Covar[c(1,4:7)]
    
    print(Covar2[s])
    
    TestCovar <- setdiff(Covar3, Covar2[s])
    
    Formula = as.formula(paste0(Resps[r], 
                                " ~ ",
                                paste(TestCovar, collapse = " + ")))
    
    TensorDevList[[Resps[r]]][[Covar2[s]]] <- bam(Formula,
                                                  data = DataList[[Resps[r]]], 
                                                  family = binomial(),
                                                  paraPen = PPList[[Resps[r]]])
    
  }
}

save(TensorDevList, file = "Output Files/TensorDevList.Rdata")
