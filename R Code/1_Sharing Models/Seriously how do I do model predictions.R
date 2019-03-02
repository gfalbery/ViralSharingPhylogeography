
# Rscript "R Code/1_Sharing Models/1b_Frequentist GAM Output.R" ####

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else 
  source("R Code/00_Master Code.R")

library(mgcv); library(tidyverse); library(ggregplot)

load("Output Files/BAMList.Rdata")

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

DataList <- PPList <- FitList <- list()

r = 1

print(Resps[r])

DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r])) %>% droplevels

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

PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), diag(nlevels(DataList[[Resps[r]]]$Sp))))


# Model Checking ####

qplot(BAMList[[Resps[r]]]$coef)

SpCoefNames <- names(BAMList[[Resps[r]]]$coef)[substr(names(BAMList[[Resps[r]]]$coef),1,5)=="SppSp"]
SpCoef <- BAMList[[Resps[r]]]$coef[SpCoefNames]
qplot(SpCoef)

# Effects ####

SpaceRange <- seq(from = min(DataList[[Resps[r]]]$Space),
                  to = max(DataList[[Resps[r]]]$Space),
                  length = 100)

PhyloRange <- #seq(from = min(scale(DataList[[Resps[r]]]$Phylo2)),
  #    to = max(scale(DataList[[Resps[r]]]$Phylo2)),
  #    length = 100)  %>% 
  0

DietRange <- #seq(from = min(DataList[[Resps[r]]]$DietSim),
  #     to = max(DataList[[Resps[r]]]$DietSim),
  #     length = 10)  %>% 
  mean(DataList[[Resps[r]]]$DietSim)

PredDF <- expand.grid(Space = SpaceRange,
                                   Phylo2 = PhyloRange,
                                   DietSim = DietRange,
                                   MinCites = mean(DataList[[Resps[r]]]$MinCites),
                                   Domestic = 0
) %>% mutate(
  
  Mean = ifelse(Space == mean(DataList[[Resps[r]]]$Space) | 
                  Phylo2 == mean(scale(DataList[[Resps[r]]]$Phylo2)), 1, 0)
  
)

PredDF$Spp <- matrix(0 , nrow = nrow(PredDF), ncol = length(SpCoef))# %>% as("dgCMatrix")

#PredDF <- PredDF %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
#                                                      PhyloQ = cut(Phylo2, quantile(Phylo2, 0:10/10),include.lowest = T, labels = 1:10))

FitPredictions  <- predict.bam(BAMList[[1]], 
                               #type = "iterms",
                               #exclude = "Spp",
                               newdata = PredDF)

FitPredictions <- FitPredictions %>% as.data.frame() %>% 
  mutate(Intercept = attr(FitPredictions,"constant"))

FitList[[Resps[r]]][,"Fit"] <- logistic(rowSums(FitPredictions))
