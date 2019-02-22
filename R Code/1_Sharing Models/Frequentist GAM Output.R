
# Frequentist GAM Output ####

SpaceRange <- seq(from = min(FinalHostMatrix$Space),
                  to = max(FinalHostMatrix$Space),
                  length = 100)

PhyloRange <- seq(from = min(FinalHostMatrix$Phylo2),
                  to = max(FinalHostMatrix$Phylo2),
                  length = 100)

DietRange <- seq(from = min(FinalHostMatrix$DietSim),
                 to = max(FinalHostMatrix$DietSim),
                 length = 100)

GAMPredDF <- expand.grid(Space = SpaceRange,
                         Phylo2 = PhyloRange,
                         DietSim = DietRange
                         #Eaten = c(0,1)
)

GAMPredDF <- GAMPredDF %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                  PhyloQ = cut(Phylo2, quantile(Phylo2, 0:10/10),include.lowest = T, labels = 1:10))

GAMPredDF[,paste0(Resps,"Fit")] <- lapply(BAMList, function(a) logistic(predict(a, newdata = GAMPredDF))) %>% bind_cols

FitSlopeTime <- gather(GAMPredDF, key = "Model", value = "Estimate", paste0(Resps[c(1,3,4,5)],"Fit"))

ggplot(FitSlopeTime %>% filter(DietSim == 0), aes(Space, Estimate)) + 
  geom_line(aes(colour = Phylo2, group = as.factor(Phylo2)), alpha = 0.5) +
  labs(y = "Predicted Value") +
  facet_wrap(~Model)

ggplot(FitSlopeTime %>% filter(DietSim == 0), aes(Phylo2, Estimate)) + 
  geom_line(aes(colour = Space, group = as.factor(Space))) +
  labs(y = "Predicted Value") +
  facet_wrap(~Model)

ggplot(FitSlopeTime %>% filter(Space == min(Space)), aes(DietSim, Estimate)) + 
  geom_line(aes(colour = Phylo2, group = as.factor(Phylo2))) +
  labs(y = "Predicted Value") +
  facet_wrap(~Model)

list(
  ggplot(GAMPredDF %>% filter(DietSim==0), aes(Phylo2, Estimate, colour = Space)) + 
    geom_line(aes(group = as.factor(Space)), alpha = 0.3),
  
  ggplot(GAMPredDF %>% filter(DietSim==0), aes(Space, Estimate, colour = Phylo2)) + 
    geom_line(aes(group = as.factor(Phylo2)), alpha = 0.3)
  
) %>% arrange_ggplot2


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
