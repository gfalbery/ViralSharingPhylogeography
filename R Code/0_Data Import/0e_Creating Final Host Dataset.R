# Creating final dataset

library(tidyverse)

rownames(Hosts) = Hosts$Sp

tCytBMatrix <- 1 - (CytBMatrix - min(CytBMatrix))/max(CytBMatrix)
tSTMatrix <- 1 - (STMatrix - min(STMatrix))/max(STMatrix)

FinalHostNames <- reduce(list(as.character(Hosts$Sp), 
                              rownames(RangeAdj1), 
                              #colnames(CytBMatrix),
                              colnames(STMatrix),
                              rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

HostThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

UpperHosts <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           PropVirus = c(HostAdj2[FHN, FHN]),
                           PropVirus2 = c(HostAdj3[FHN, FHN]),
                           Space = c(RangeAdj1[FHN, FHN]),
                           #Phylo = c(tCytBMatrix[FHN, FHN]), 
                           Phylo2 = c(tSTMatrix[FHN, FHN]),
                           Sp = as.character(rep(FHN, each = length(FHN))),
                           Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

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
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp2 + 1)
HostMatrixdf$MinCites <- apply(HostMatrixdf[,c("hAllZACites", "hAllZACites.Sp2")],1, function(a) min(a, na.rm = T))

HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$MinDCites <- apply(HostMatrixdf[,c("hDiseaseZACites", "hDiseaseZACites.Sp2")],1, function(a) min(a, na.rm = T))
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp2 + 1)

HostMatrixdf$SpaceQuantile <- cut(HostMatrixdf$Space, 
                                  breaks = c(-0.1,0,0.25,0.5, 0.75, 1.1),
                                  labels = c(0, 0.25, 0.5, 0.75, 1))

#HostMatrixdf$PhyloQuantile <- cut(HostMatrixdf$Phylo2, 
#                                  breaks = c(quantile(HostMatrixdf$Phylo2, 0:5/5)),
#                                  labels = c(0, 0.25, 0.5, 0.75, 1))

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostThemselves <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo2
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)
FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))

# Establishing all mammal data ####

# Replacing absent names in the full ST matrix ####

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]
AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

NameReplace <- c(
  "Micaelamys_namaquensis",
  "Akodon_paranaensis",
  "Bos_frontalis",
  "Bos_grunniens",
  "Bubalus_arnee", # Absent
  "Capra_hircus",
  "Hexaprotodon_liberiensis",
  "Equus_burchellii",
  "Oryzomys_alfaroi" ,
  "Oryzomys_laticeps",
  "Oryzomys_megacephalus",
  "Callithrix_argentata",
  "Miniopterus_schreibersii",
  "Myotis_ricketti",
  "Oryzomys_albigularis",
  "Ovis_aries",
  "Piliocolobus_badius",
  "Piliocolobus_rufomitratus" ,
  "Lycalopex_gymnocercus" ,
  "Rhinolophus_hildebrandtii",
  "Oryzomys_angouya",
  "Mops_condylurus",
  "Chaerephon_plicatus",
  "Chaerephon_pumilus",
  "Taurotragus_oryx")

names(NameReplace) <- AbsentHosts

rownames(FullSTMatrix) <- colnames(FullSTMatrix) <- sapply(rownames(FullSTMatrix), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

