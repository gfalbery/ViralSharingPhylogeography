

HostThems <- which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag = T))

MammalThems <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T)&
                       lower.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))


ObsFullCors <- HostMatrixdf %>% slice(-HostThems) %>%
  group_by(Sp) %>% summarise(SpaceMean = mean(Space),
                             PhyloMean = mean(Phylo),
                             SpaceSum = sum(Space),
                             PhyloSum = sum(Phylo)) %>%
  left_join(
    
    AllMammalMatrix %>% slice(-MammalThems) %>%
      group_by(Sp) %>% summarise(SpaceMean = mean(Space),
                                 PhyloMean = mean(Phylo),
                                 SpaceSum = sum(Space),
                                 PhyloSum = sum(Phylo)),
    
    by = "Sp", suffix = c(".Obs", ".Full")
    
  ) 

jpeg("SIFigures/ObsKnownCorrelations.jpeg", units = "mm", width = 180, height = 100, res = 300)

list(
ggplot(ObsFullCors, aes(SpaceMean.Obs, SpaceMean.Full)) + 
  geom_point() +
  labs(x = "Observed Mean Space Sharing",
       y = "Global Mean Space Sharing") +
  geom_smooth(colour = "black"),

ggplot(ObsFullCors, aes(PhyloMean.Obs, PhyloMean.Full)) + 
  geom_point() +
  labs(x = "Observed Mean Phylogenetic Relatedness",
       y = "Global Mean Phylogenetic Relatedness") +
  geom_smooth(colour = "black")

) %>% arrange_ggplot2

dev.off()

ggplot(ObsFullCors, aes(SpaceSum.Obs, SpaceSum.Full)) + 
  geom_point() +
  geom_smooth()

ggplot(ObsFullCors, aes(PhyloSum.Obs, PhyloSum.Full)) + 
  geom_point() +
  geom_smooth()
