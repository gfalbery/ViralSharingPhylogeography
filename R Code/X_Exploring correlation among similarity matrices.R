# X_ Exploring Similarity Matrix Models

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

library(GGally)

HostThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

png(filename = "Figures/Pairwise Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
ggpairs(LongMatrixdf[-HostThemselves,], 
                lower = list(continuous = wrap("smooth", method = "gam"))) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space and Phylogeny Correlate with Viral Sharing")
dev.off()

# For Viruses ####

VirusThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(VirusAdj[FVN,FVN], diag = T) &
          lower.tri(VirusAdj[FVN,FVN], diag  = T))

png(filename = "Figures/Pairwise Virus Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
ggpairs(VirusLongMatrixdf[-VirusThemselves,], 
                lower = list(continuous = wrap("smooth", method = "gam"))) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space Correlates with Host Sharing")
dev.off()

