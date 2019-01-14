# Creating final dataset

library(tidyverse)

FinalHostNames <- reduce(list(as.character(SpatialHosts$Sp), 
                              rownames(RangeAdj1), 
                              colnames(CytBMatrix),
                              rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

LongMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           PropVirus = c(HostAdj2[FHN, FHN]),
                           Space = c(RangeAdj1[FHN, FHN]),
                           Phylo = c(1-CytBMatrix[FHN, FHN]) # Gonna invert this
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

IM1 <- inla(data = LongMatrixdf, 
            PropVirus ~ Space + Phylo,
            family = "gamma")

summary(IM1)

Efxplot(list(IM1))

IM1 <- inla(data = LongMatrixdf, # Doesn't fit
            Virus ~ Space + Phylo,
            control.compute = list(dic = TRUE),
            family = "nbinomial")

IM2 <- inla(data = LongMatrixdf, # Doesn't fit
            Virus ~ Space + Phylo,
            control.compute = list(dic = TRUE),
            family = "zeroinflatednbinomial1")

mf = 1

MC1 <- MCMCglmm(data = LongMatrixdf, 
                Virus ~ Space + Phylo,
                family = "poisson",
                nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
                thin = 10*mf,burnin=3000*mf)

MC2 <- MCMCglmm(data = LongMatrixdf, 
                Virus ~ Space + Phylo,
                rcov =~ idh(trait):units, 
                family = "zipoisson",
                nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
                thin = 10*mf,burnin=3000*mf)

MC3 <- MCMCglmm(data = LongMatrixdf, 
                Virus ~ trait -1 + trait:(Space + Phylo),
                rcov =~ idh(trait):units, 
                family = "zipoisson",
                nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
                thin = 10*mf,burnin=3000*mf)

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
  which(upper.tri(VirusAdj[FVN,FVN], diag = T)&
          lower.tri(VirusAdj[FVN,FVN], diag  = T))

png(filename = "Figures/Virus Pairwise Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
GGally::ggpairs(VirusLongMatrixdf[-VirusThemselves,], 
                lower = list(continuous = "smooth", method = "gam")) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space Correlates with Host Sharing")
dev.off()

FVN2 <- reduce(list(SpatialViruses$Sp, 
                    colnames(VirusHostPD)[-which(is.na(VirusHostPD[1,]))],
                    rownames(VirusRangeAdj1)[which(sapply(GridList[unique(SpatialViruses$Sp)], length)>0)]), 
               intersect)



