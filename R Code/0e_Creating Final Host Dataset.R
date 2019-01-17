# Creating final dataset

library(tidyverse)

tCytBMatrix <- 1 - (CytBMatrix - min(CytBMatrix))/max(CytBMatrix)

FinalHostNames <- reduce(list(as.character(Hosts$Sp), 
                              rownames(RangeAdj1), 
                              colnames(CytBMatrix),
                              rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

HostMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           PropVirus = c(HostAdj2[FHN, FHN]),
                           PropVirus2 = c(HostAdj3[FHN, FHN]),
                           Space = c(RangeAdj1[FHN, FHN]),
                           Phylo = c(tCytBMatrix[FHN, FHN]), # Gonna invert this
                           Host = rep(FHN, each = length(FHN)) 
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



