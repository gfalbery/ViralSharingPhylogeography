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
                           Sp = rep(FHN, each = length(FHN)),
                           Sp2 = rep(FHN, length(FHN))
)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites", 
                   "LongMean", "LatMean")

HostMatrixdf[,HostMatrixVar] <- Hosts[HostMatrixdf$Sp, HostMatrixVar]
HostMatrixdf[,paste0(HostMatrixVar,".Sp2")] <- Hosts[HostMatrixdf$Sp2, HostMatrixVar]
HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] <- HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] %>% mutate(hAllZACites = 1167, hDiseaseZACites = 115)

HostMatrixdf <- HostMatrixdf %>% mutate(
  hOrder = Hosts[HostMatrixdf$Sp,"hOrder"],
  hFamily = Hosts[HostMatrixdf$Sp,"hFamily"],
  hDom = Hosts[HostMatrixdf$Sp,"hDom"]
)

HostMatrixdf$Space0 <- ifelse(HostMatrixdf$Space == 0, "No Overlap", "Overlap")
HostMatrixdf$Cites <- log(HostMatrixdf$hAllZACites + 1)
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp + 1)
HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp + 1)

HostMatrixdf$SpaceQuantile <- cut(HostMatrixdf$Space, 
                                  breaks = c(-0.1,0,0.25,0.5, 0.75, 1.1),
                                  labels = c(0, 0.25, 0.5, 0.75, 1))

HostMatrixdf$PhyloQuantile <- cut(HostMatrixdf$Phylo, 
                                  breaks = c(-0.1,0,0.25,0.5, 0.75, 1.1),
                                  labels = c(0, 0.25, 0.5, 0.75, 1))

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



