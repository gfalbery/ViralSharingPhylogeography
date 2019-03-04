
# Rscript "R Code/1_Sharing Models/1b_Frequentist GAM Output.R" ####

if(file.exists("Output Files/Finaldf.Rdata")) load("Output Files/Finaldf.Rdata") else source("R Code/00_Master Code.R")

library(mgcv); library(tidyverse); library(ggregplot)

load("Output Files/BAMList.Rdata")
load("Output Files/BAMList2.Rdata")

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

DataList <- PPList <- FitList <- list()

for(r in 1:length(BAMList)){
  
  print(Resps[r])
  
  # Model Checking ####
  
  qplot(BAMList[[Resps[r]]]$coef)
  
  SpCoefNames <- names(BAMList[[Resps[r]]]$coef)[substr(names(BAMList[[Resps[r]]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[Resps[r]]]$coef[SpCoefNames]

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
  
  FitList[[Resps[r]]] <- expand.grid(Space = SpaceRange,
                                     Phylo = PhyloRange,
                                     DietSim = DietRange,
                                     MinCites = mean(DataList[[Resps[r]]]$MinCites),
                                     Domestic = 0
  )
  
  FitList[[Resps[r]]]$Spp <- matrix(0 , nrow = nrow(FitList[[Resps[r]]]), ncol = length(SpCoef))# %>% as("dgCMatrix")
  
  FitList[[Resps[r]]] <- FitList[[Resps[r]]] %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                                        PhyloQ = cut(Phylo, quantile(Phylo, 0:10/10),include.lowest = T, labels = 1:10))
  
  FitPredictions  <- predict.gam(BAMList[[1]], 
                                 newdata = FitList[[Resps[r]]])
  
  #FitPredictions <- FitPredictions %>% as.data.frame() %>% 
  #  mutate(Intercept = attr(FitPredictions,"constant"))
  
  #FitList[[Resps[r]]][,"Fit"] <- FitPredictions
  FitList[[Resps[r]]][,"Fit"] <- logistic(FitPredictions)
  
}

save(FitList, file = "Output Files/FitList.Rdata")

r=1

FitList[["VirusBinary"]] %>% filter(Phylo == last(PhyloRange)) %>%
  ggplot(aes(Space, Fit)) + geom_line() +
  lims(y = c(0,1))

FitList[["VirusBinary"]] %>% filter(Space == SpaceRange[10]) %>%
  ggplot(aes(Phylo, Fit)) + geom_line() +
  lims(y = c(0,1))

FitList[["VirusBinary"]] %>% filter(Space == last(SpaceRange)) %>%
  ggplot(aes(Phylo, Fit)) + geom_line()

FitList[["VirusBinary"]] %>% filter(Phylo == last(PhyloRange)) %>%
  ggplot(aes(Space, Fit)) + geom_line()

tiff("Figures/Model Predictions.jpeg", units = "mm", width = 200, height = 150, res = 300)

list(
  ggplot(FitList[["VirusBinary"]], aes(Phylo, Fit, colour = Space)) + 
    geom_line(aes(group = as.factor(Space)), alpha = 0.3) +
    geom_rug(data = DataList[[Resps[r]]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
  
  ggplot(FitList[["VirusBinary"]], aes(Space, Fit, colour = Phylo)) + 
    geom_line(aes(group = as.factor(Phylo)), alpha = 0.3) +
    geom_rug(data = DataList[[Resps[r]]], inherit.aes = F, aes(x = Space), alpha = 0.01)
  
) %>% arrange_ggplot2

dev.off()

