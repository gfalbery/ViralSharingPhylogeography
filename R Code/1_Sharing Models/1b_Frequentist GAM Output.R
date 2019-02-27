
# Frequentist GAM Output ####

library(mgcv); library(tidyverse); library(ggregplot)

load("Output Files/BAMList.Rdata")

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

DataList <- PPList <- list()

for(r in 1:length(BAMList)){
  
  print(Resps[r])
  
  DataList[[Resps[r]]] <- FinalHostMatrix %>% filter(!is.na(Resps[r]))
  
  DataList[[Resps[r]]]$Sp <- factor(DataList[[Resps[r]]]$Sp, levels = union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2))
  DataList[[Resps[r]]]$Sp2 <- factor(DataList[[Resps[r]]]$Sp2, levels = union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2))
  
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

# Model Checking ####

qplot(BAMList[[1]]$coef)

SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
SpCoef <- BAMList[[1]]$coef[SpCoefNames]
qplot(SpCoef)

# Effects ####

SpaceRange <- seq(from = min(FinalHostMatrix$Space),
                  to = max(FinalHostMatrix$Space),
                  length = 20) %>% c(mean(FinalHostMatrix$Space))

PhyloRange <- seq(from = min(scale(FinalHostMatrix$Phylo2)),
                  to = max(scale(FinalHostMatrix$Phylo2)),
                  length = 20)  %>% c(mean(scale(FinalHostMatrix$Phylo2)))

DietRange <- seq(from = min(FinalHostMatrix$DietSim),
                 to = max(FinalHostMatrix$DietSim),
                 length = 10)  %>% c(mean(FinalHostMatrix$DietSim))

GAMPredDF <- expand.grid(Space = SpaceRange,
                         Phylo2 = PhyloRange,
                         DietSim = DietRange,
                         MinCites = mean(FinalHostMatrix$MinCites),
                         Domestic = 0
) %>% mutate(
  
  Mean = ifelse(Space == mean(FinalHostMatrix$Space) | 
                  Phylo2 == mean(scale(FinalHostMatrix$Phylo2)), 1, 0)
  
)

GAMPredDF$Spp <- matrix(0 , nrow = nrow(GAMPredDF), ncol = length(SpCoef))# %>% as("dgCMatrix")

GAMPredDF <- GAMPredDF %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                  PhyloQ = cut(Phylo2, quantile(Phylo2, 0:10/10),include.lowest = T, labels = 1:10))

FitPredictions  <- predict.gam(BAMList[[1]], 
                               type = "terms",
                               newdata = GAMPredDF)

FitPredictions <- FitPredictions %>% as.data.frame() %>% 
  mutate(Intercept = attr(FitPredictions,"constant"))

GAMPredDF[,"Fit"] <- logistic(rowSums(FitPredictions))

tiff("Figures/Model Predictions.jpeg", units = "mm", width = 200, height = 150, res = 300)

list(
  ggplot(GAMPredDF %>% filter(DietSim==0), aes(Phylo2, Fit, colour = Space)) + 
    geom_line(aes(group = as.factor(Space), lty = as.factor(Mean)), alpha = 0.3),
  
  ggplot(GAMPredDF %>% filter(DietSim==0), aes(Space, Fit, colour = Phylo2)) + 
    geom_line(aes(group = as.factor(Phylo2), lty = as.factor(Mean)), alpha = 0.3)
  
) %>% arrange_ggplot2

dev.off()

list(
  ggplot(GAMPredDF %>% filter(DietSim==0), aes(Phylo2, Estimate, colour = Space)) + 
    geom_line(aes(group = as.factor(Space)), alpha = 0.3) +
    geom_rug(data = FinalHostMatrix, inherit.aes = F, aes(x = Phylo2), alpha = 0.01),
  
  ggplot(GAMPredDF %>% filter(DietSim==0), aes(Space, Estimate, colour = Phylo2)) + 
    geom_line(aes(group = as.factor(Phylo2)), alpha = 0.3) +
    geom_rug(data = FinalHostMatrix, inherit.aes = F, aes(x = Space), alpha = 0.01)
  
) %>% arrange_ggplot2


BinBAM2 <- bam(VirusBinary ~ t2(Space, scale(Phylo2)) + t2(Space, DietSim),# + s(DietSim),# + Eaten, # + 
               #(1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
               data = FinalHostMatrix, 
               family = binomial())

summary(BinBAM2)
