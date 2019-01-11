# Creating final dataset

library(tidyverse)

FinalHostNames <- reduce(list(SpatialHosts$Sp, 
                              rownames(RangeAdj1), 
                              rownames(CytBMatrix),
                              rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

LongMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           PropVirus = c(HostAdj2[FHN, FHN]),
                           Space = c(RangeAdj1[FHN, FHN]),
                           Phylo = 1-c(CytBMatrix[FHN, FHN]) # Gonna invert this
)

Themselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

png(filename = "Figures/Pairwise Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
GGally::ggpairs(LongMatrixdf[-Themselves,], 
                lower = list(continuous = "smooth", method = "gam")) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space and Phylogeny Correlate with Viral Sharing")
dev.off()

# Makes sense (NB Phylo has been inverted so no longer a measure of distance).
# Better correlations with absolute virus count than observation-weighted count.

FHosts <- SpatialHosts[SpatialHosts$Sp%in%FHN,]
FHosts <- FHosts[order(FHosts$Sp),]
FHosts$Pixels <- diag(RangeOverlap)[FHN]
FHosts$ViralRichness <- diag(HostAdj)[FHN]

# Virus dataset ####


FinalVirusNames <- reduce(list(SpatialViruses$Sp, 
                               rownames(VirusRangeAdj1)[which(sapply(GridList[unique(SpatialViruses$Sp)], length)>0)]), 
                          intersect)

FVN <- FinalVirusNames; length(FVN)

FViruses <- SpatialViruses[SpatialViruses$Sp%in%FVN,]
FViruses <- FViruses[order(FViruses$Sp),]
FViruses$Pixels <- diag(VirusRangeOverlap)[FVN]
FViruses$ViralRichness <- diag(VirusRangeAdj1)[FVN]


VirusLongMatrixdf <- data.frame(Host = c(VirusAdj[FVN, FVN]),
                           PropHost = c(VirusAdj2[FVN, FVN]),
                           Space = c(VirusRangeAdj1[FVN, FVN]) # Gonna invert this
)

VirusThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(VirusAdj[FVN,FVN], diag = T)&lower.tri(VirusAdj[FVN,FVN], diag  = T))

png(filename = "Figures/Virus Pairwise Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
GGally::ggpairs(VirusLongMatrixdf[-VirusThemselves,], 
                lower = list(continuous = "smooth", method = "gam")) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space Does Not Correlate with Host Sharing")
dev.off()




