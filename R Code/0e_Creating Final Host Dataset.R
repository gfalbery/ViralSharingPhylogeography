# Creating final dataset

library(tidyverse)

rownames(Hosts) = Hosts$Sp

tCytBMatrix <- 1 - (CytBMatrix - min(CytBMatrix))/max(CytBMatrix)

FinalHostNames <- reduce(list(as.character(Hosts$Sp), 
                              rownames(RangeAdj1), 
                              colnames(CytBMatrix),
                              rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

HostThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

HostMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           PropVirus = c(HostAdj2[FHN, FHN]),
                           PropVirus2 = c(HostAdj3[FHN, FHN]),
                           Space = c(RangeAdj1[FHN, FHN]),
                           Phylo = c(tCytBMatrix[FHN, FHN]), # Gonna invert this
                           Sp = rep(FHN, each = length(FHN)) 
)

HostMatrixdf <- HostMatrixdf %>% mutate(
  hOrder = Hosts[HostMatrixdf$Sp,"hOrder"],
  hFamily = Hosts[HostMatrixdf$Sp,"hFamily"],
  hDom = Hosts[HostMatrixdf$Sp,"hDom"]
)

FHosts <- Hosts[Hosts$Sp%in%FHN,]
FHosts <- FHosts[order(FHosts$Sp),]

# Virus dataset ####

FinalVirusNames <- reduce(list(Viruses$Sp, 
                               rownames(VirusRangeAdj1)[which(sapply(GridList[unique(Viruses$Sp)], length)>0)]), 
                          intersect)

FVN <- FinalVirusNames; length(FVN)

FViruses <- Viruses[Viruses$Sp%in%FVN,]
FViruses <- FViruses[order(FViruses$Sp),]
FViruses$Pixels <- diag(VirusRangeOverlap)[FVN]
FViruses$ViralRichness <- diag(VirusRangeAdj1)[FVN]

VirusLongMatrixdf <- data.frame(Host = c(VirusAdj[FVN, FVN]),
                                PropHost = c(VirusAdj2[FVN, FVN]),
                                Space = c(VirusRangeAdj1[FVN, FVN]) # Gonna invert this
)

# Plotting some exploratory stuff ####

ggplot(HostMatrixdf[-HostThemselves,], aes(Space, Virus)) + 
  geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(aes(colour = hDom)) +
  stat_smooth(geom = "ribbon", fill = NA, lty = 2, aes(group = hDom))

ggplot(HostMatrixdf[-HostThemselves,], aes(Phylo, Virus)) + 
  geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(aes(colour = hDom)) +
  stat_smooth(geom = "ribbon", fill = NA, lty = 2, aes(group = hDom))

ggplot(HostMatrixdf[-HostThemselves,], aes(Space, Virus)) + 
  geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(aes(colour = hOrder)) +
  stat_smooth(geom = "ribbon", fill = NA, lty = 2, aes(group = hOrder))

ggplot(HostMatrixdf[-HostThemselves,], aes(Phylo, Virus)) + 
  geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(aes(colour = hDom)) +
  stat_smooth(geom = "ribbon", fill = NA, lty = 2, aes(group = hDom))



