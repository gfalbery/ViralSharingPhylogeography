# Creating final dataset

FinalHostNames <- intersect(Hosts$Sp, rownames(RangeAdj1))
FinalHostNames <- intersect(FinalHostNames, rownames(CytBMatrix))
FinalHostNames <- intersect(FinalHostNames, rownames(HostAdj))

FHN <- FinalHostNames; length(FHN)

LongMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           PropVirus = c(HostAdj2[FHN, FHN]),
                           Space = c(RangeAdj1[FHN, FHN]),
                           Phylo = -c(CytBMatrix[FHN, FHN]) # Gonna invert this
)

Themselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

GGally::ggpairs(LongMatrixdf[-Themselves,], 
                lower = list(continuous = "smooth"),
                method = "gam") + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey"))

# Makes sense (NB Phylo has been converted from a measure of distance).
# Better correlations with absolute virus count than observation-weighted count.

FHosts <- Hosts[Hosts$Sp%in%FHN,]
FHosts <- FHosts[order(FHosts$Sp),]
FHosts$Pixels <- diag(RangeOverlap)[FHN]
FHosts$ViralRichness <- diag(HostAdj)[FHN]


