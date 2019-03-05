# Extracting deviance parameters from BAM ####

# Rscript "R Code/1_Sharing Models/1c_Frequentist GAM Deviance Contributions.R" 

# Running Frequentist GAMS

#if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

load("Output Files/BAMList.Rdata")
#load("Output Files/BAMList2.Rdata")

library(mgcv); library(tidyverse)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

DevList <- DevList2 <- DataList <- PPList <- list()

for(r in 1:length(BAMList)){
  
  print(Resps[r])
  
  Covar <- c(#"s(Space)", 
    #"s(scale(Phylo))",
    "t2(Space,Phylo)",
    "s(DietSim)",
    "MinCites", 
    "Domestic",
    "Spp")
  
  Formula = as.formula(paste0(Resps[r], 
                              " ~ ",
                              paste(Covar, collapse = " + ")))
  
  DevList[[Resps[r]]]$FullModel <- bam(Formula,
                                       data = DataList[[Resps[1]]], 
                                       family = binomial(),
                                       paraPen = PPList[[Resps[1]]])
  
  for(s in 1:length(Covar)){
    
    print(Covar[s])
    
    TestCovar <- setdiff(Covar, Covar[s])
    
    Formula = as.formula(paste0(Resps[r], 
                                " ~ ",
                                paste(TestCovar, collapse = " + ")))
    
    DevList[[Resps[r]]][[Covar[s]]] <- bam(Formula,
                                           data = DataList[[Resps[r]]], 
                                           family = binomial(),
                                           paraPen = PPList[[Resps[r]]])
    
  }
}

save(DevList, file = "Output Files/DevList.Rdata")
