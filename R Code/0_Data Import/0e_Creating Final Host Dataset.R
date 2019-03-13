# Creating final dataset

# Replacing absent names in the full ST matrix ####

FinalHostNames <- reduce(list(
  rownames(RangeAdj), 
  colnames(STMatrix),
  rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

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

#rownames(VD) <- colnames(VD) <- sapply(rownames(VD), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

# Going Ahead ####

library(tidyverse)

rownames(Hosts) = Hosts$Sp

tCytBMatrix <- 1 - (CytBMatrix - min(CytBMatrix))/max(CytBMatrix)
tSTMatrix <- 1 - (STMatrix - min(STMatrix))/max(STMatrix)

FinalHostNames <- reduce(list(
  rownames(RangeAdj), 
  colnames(STMatrix),
  #colnames(VD),
  rownames(HostAdj)), intersect)

FHN <- FinalHostNames; length(FHN)

UpperHosts <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           Space = c(RangeAdj[FHN, FHN]),
                           Phylo = c(tSTMatrix[FHN, FHN]),
                           Sp = as.character(rep(FHN, each = length(FHN))),
                           Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites"
                   #"LongMean", "LatMean")
)

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

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

#FinalHostMatrix <- FinalHostMatrix %>% left_join(LongDiet, by = c("Sp","Sp2")) %>%
#  mutate(DietSim = 1 - DietSim)

#FinalHostMatrix <- FinalHostMatrix %>% left_join(EltonTraits[,c("Scientific","Carnivore")], by = c("Sp" = "Scientific"))
#FinalHostMatrix <- FinalHostMatrix %>% left_join(EltonTraits[,c("Scientific","Carnivore")], by = c("Sp2" = "Scientific"),
#                                                 suffix = c("",".Sp2"))

#FinalHostMatrix$Eaten <- ifelse(FinalHostMatrix$Carnivore==FinalHostMatrix$Carnivore.Sp2,0,1)

Remove1 <- FinalHostMatrix %>% group_by(Sp) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% select(Sp)
Remove2 <- FinalHostMatrix %>% group_by(Sp2) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% select(Sp2)

Remove3 <- which(table(c((FinalHostMatrix %>% filter(Phylo < 0.25) %>% select(Sp, Sp2))$Sp %>% as.character(),
                         (FinalHostMatrix %>% filter(Phylo < 0.25) %>% select(Sp, Sp2))$Sp2 %>% as.character()))>20) %>% 
  names

RemoveSp <- union(intersect(Remove1$Sp, Remove2$Sp2), Remove3)

FinalHostMatrix <- FinalHostMatrix %>% filter(!Sp%in%RemoveSp&!Sp2%in%RemoveSp)

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

FinalHostMatrix <- FinalHostMatrix %>% slice(order(Sp,Sp2))

# Let's try this a second time ####

FHN <- levels(FinalHostMatrix$Sp)

HostMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           Space = c(RangeAdj[FHN, FHN]),
                           Phylo = c(tSTMatrix[FHN, FHN]),
                           Sp = as.character(rep(FHN, each = length(FHN))),
                           Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites"
                   #"LongMean", "LatMean")
)

HostMatrixdf[,HostMatrixVar] <- Hosts[HostMatrixdf$Sp, HostMatrixVar]
HostMatrixdf[,paste0(HostMatrixVar,".Sp2")] <- Hosts[HostMatrixdf$Sp2, HostMatrixVar]

HostMatrixdf$Cites <- log(HostMatrixdf$hAllZACites + 1)
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp2 + 1)
HostMatrixdf$MinCites <- apply(HostMatrixdf[,c("hAllZACites", "hAllZACites.Sp2")],1, function(a) min(a, na.rm = T))

HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$MinDCites <- apply(HostMatrixdf[,c("hDiseaseZACites", "hDiseaseZACites.Sp2")],1, function(a) min(a, na.rm = T))
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp2 + 1)

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

FinalHostMatrix$Gz <- as.numeric(FinalHostMatrix$Space==0)

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

FinalHostMatrix <- FinalHostMatrix %>% slice(order(Sp,Sp2)) %>%
  mutate(Phylo = (Phylo - min(Phylo))/max(Phylo - min(Phylo)))
