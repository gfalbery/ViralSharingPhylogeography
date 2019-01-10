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

png(filename = "Figures/Pairwise.jpg", units = "mm", width = 200, height = 200, res = 300)
GGally::ggpairs(LongMatrixdf[-Themselves,], 
                lower = list(continuous = "smooth", method = "gam")) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space and Phylogeny Correlate with Viral Sharing")
dev.off()

# Makes sense (NB Phylo has been inverted so no longer a measure of distance).
# Better correlations with absolute virus count than observation-weighted count.

FHosts <- Hosts[Hosts$Sp%in%FHN,]
FHosts <- FHosts[order(FHosts$Sp),]
FHosts$Pixels <- diag(RangeOverlap)[FHN]
FHosts$ViralRichness <- diag(HostAdj)[FHN]


