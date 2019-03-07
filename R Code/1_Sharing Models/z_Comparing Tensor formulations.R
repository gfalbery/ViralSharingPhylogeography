
# Trying different tensor combos

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else{
  source("R Code/00_Master Code.R")
}

library(mgcv); library(tidyverse); library(ggregplot)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")[1]

TryList <- DataList <- PPList <- FormulaList <- list()

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
  
  PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), 
                                        diag(nlevels(DataList[[Resps[r]]]$Sp))))
  
  Covar <- c("t2(Space, Phylo)",
             "ti(Space, Phylo)",
             "s(Space)",
             "s(Phylo)",
             "s(DietSim)",
             "MinCites",
             "Domestic",
             "Spp")
  
  print("T2")
  
  FormulaList$T2 = as.formula(paste0(Resps[r], 
                                     " ~ ",
                                     paste(Covar[c(1,5:8)], collapse = " + ")
  ))
  
  # Tensor
  
  TryList[[Resps[r]]]$T2 <- bam(FormulaList$T2,
                                data = DataList[[Resps[r]]], 
                                family = binomial(),
                                paraPen = PPList[[Resps[r]]], select = T
  )
  
  FormulaList$Te = as.formula(paste0(Resps[r], 
                                     " ~ ",
                                     paste(Covar[c(2:8)], collapse = " + ")
  ))
  
  print("Te1")
  
  TryList[[Resps[r]]]$Te1 <- bam(FormulaList$Te,
                                 data = DataList[[Resps[r]]], 
                                 family = binomial(),
                                 paraPen = PPList[[Resps[r]]]
  )
  
  print("Te2")
  
  TryList[[Resps[r]]]$Te2 <- bam(FormulaList$Te,
                                 data = DataList[[Resps[r]]], 
                                 family = binomial(),
                                 paraPen = PPList[[Resps[r]]], select = T
  )
  
  KnotCovar <- c(#"t2(Space, Phylo)",
    "ti(Space, Phylo, k = 20)",
    "s(Space, k = 12)",
    "s(Phylo, k = 12)")
  
  for(Knots in 1:length(KnotCovar)){
    
    print(KnotCovar[Knots])
    
    FormulaList[[KnotCovar[Knots]]] = as.formula(paste0(Resps[r], 
                                                        " ~ ",
                                                        paste(Covar[c(2:8)][-Knots], collapse = " + "),
                                                        " + ", KnotCovar[Knots]
    ))
    
    TryList[[Resps[r]]][[KnotCovar[Knots]]] <- bam(FormulaList[[KnotCovar[Knots]]],
                                                   data = DataList[[Resps[r]]], 
                                                   family = binomial(),
                                                   paraPen = PPList[[Resps[r]]], select = T
    )
    
  }
  
  
}

save(TryList, file = "Output Files/TryList.Rdata")

sapply(TryList[[Resps[r]]], AIC)

TryFitList <- list()

Tries <- c("T2", "Te1", "Te2", KnotCovar)

for(s in 1:length(TryList)){
  
  Model <- TryList[[Resps[r]]][[Tries[s]]]
  
  print(Tries[s])
  
  # Model Checking ####
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  # Effects ####
  
  SpaceRange <- seq(from = min(DataList[[Resps[r]]]$Space),
                    to = max(DataList[[Resps[r]]]$Space),
                    length = 100) %>% 
    c(mean(DataList[[Resps[r]]]$Space))
  
  PhyloRange <- seq(from = min(DataList[[Resps[r]]]$Phylo),
                    to = max(DataList[[Resps[r]]]$Phylo),
                    length = 100)  %>% 
    c(mean(DataList[[Resps[r]]]$Phylo))
  
  DietRange <- seq(from = min(DataList[[Resps[r]]]$DietSim),
                   to = max(DataList[[Resps[r]]]$DietSim),
                   length = 10)  %>% 
    c(mean(DataList[[Resps[r]]]$DietSim))
  
  TryFitList[[Resps[r]]][[Tries[s]]] <- expand.grid(Space = SpaceRange,
                                     Phylo = PhyloRange,
                                     DietSim = DietRange,
                                     MinCites = mean(DataList[[Resps[r]]]$MinCites),
                                     Domestic = 0
  )
  
  TryFitList[[Resps[r]]][[Tries[s]]]$Spp <- matrix(0 , nrow = nrow(TryFitList[[Resps[r]]][[Tries[s]]]), ncol = length(SpCoef))# %>% as("dgCMatrix")
  
  TryFitList[[Resps[r]]][[Tries[s]]] <- TryFitList[[Resps[r]]][[Tries[s]]] %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                                        PhyloQ = cut(Phylo, quantile(Phylo, 0:10/10),include.lowest = T, labels = 1:10))
  
  FitPredictions  <- predict.gam(Model, 
                                 newdata = TryFitList[[Resps[r]]][[Tries[s]]])
  
  TryFitList[[Resps[r]]][[Tries[s]]][,"Fit"] <- logistic(FitPredictions)
}

save(TryFitList, file = "Output Files/TryFitList.Rdata")

