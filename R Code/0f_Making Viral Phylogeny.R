library(geiger) 
library(ape)
library(picante)

vt <- read.nexus("data/ViralTree.tre")

vt$tip.label[sort.list(vt$tip.label)]  #cladogram tree, sort taxon list
TipDrop <- setdiff(vt$tip.label, Viruses$Sp) #compare two lists, elements from 1st not in 2nd. Can use to then create list and drop.tips

setdiff(Viruses$Sp, vt$tip.label)

VirusTree <- drop.tip(vt, which(vt$tip.label %in% TipDrop)) #new tree, drop all tips in cytb not in h

VirusPD <- as.data.frame(cophenetic(VirusTree)) %>% as.matrix
VirusPD <- VirusPD/(10^-314)

FinalVirusNames <- reduce(list(Viruses$Sp, 
                               rownames(VirusRangeAdj1)[which(sapply(GridList[unique(Viruses$Sp)], length)>0)],
                               rownames(VirusPD)), 
                          intersect)

FVN <- FinalVirusNames; length(FVN)

tVirusPD <- 1 - (VirusPD - min(VirusPD))/max(VirusPD)

VirusLongMatrixdf <- data.frame(Host = c(VirusAdj[FVN, FVN]),
                                PropHost = c(VirusAdj2[FVN, FVN]),
                                Space = c(VirusRangeAdj1[FVN, FVN]),
                                Phylo = c(tVirusPD[FVN,FVN])
)


VirusThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(VirusAdj[FVN,FVN], diag = T)&lower.tri(VirusAdj[FVN,FVN], diag  = T))

UpperViruses <- # Removing diagonals, as they're uninformative
  which(upper.tri(VirusAdj[FVN,FVN], diag = T))

ggpairs(VirusLongMatrixdf[-UpperViruses,], lower = list(continuous = "smooth"))
