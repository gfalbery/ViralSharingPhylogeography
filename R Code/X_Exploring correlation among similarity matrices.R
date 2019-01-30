# X_ Exploring Host Similarity Matrix Models

library(MCMCglmm)

HostMatrixdf$PropVirus[HostMatrixdf$PropVirus==0] <- 0.0001
HostMatrixdf$PropVirus[HostMatrixdf$PropVirus==1] <- 0.9999

IM1 <- inla(data = HostMatrixdf[-HostThemselves,], 
            PropVirus ~ Space + Phylo,
            family = "beta")

summary(IM1) # Interesting space is more important here ####

Efxplot(list(IM1))

IM1 <- inla(data = FinalHostMatrix, # Doesn't fit
            Virus ~ Space + Phylo2 + Space:Phylo2 + MinCites + DomDom,
            control.compute = list(dic = TRUE),
            family = "zeroinflatednbinomial1")

IM2 <- inla(data = HostMatrixdf[-HostThemselves,], # Doesn't fit
            Virus ~ Space + Phylo,
            control.compute = list(dic = TRUE),
            family = "zeroinflatednbinomial1")

mf = 1

MC1 <- MCMCglmm(data = HostMatrixdf, 
                Virus ~ Space + Phylo,
                family = "poisson",
                nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
                thin = 10*mf,burnin=3000*mf)

MC2 <- MCMCglmm(data = HostMatrixdf[-HostThemselves,], 
                Virus ~ Space + Phylo,
                rcov =~ idh(trait):units, 
                family = "zipoisson",
                nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
                thin = 10*mf, burnin=3000*mf)

MC3 <- MCMCglmm(data = HostMatrixdf[-HostThemselves,], 
                Virus ~ trait -1 + trait:(Space + Phylo),
                rcov =~ idh(trait):units, 
                family = "zipoisson",
                nitt = 13000*mf, # REMEMBER YOU'VE DONE THIS
                thin = 10*mf, burnin=3000*mf)

arrange_ggplot2(list(
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, Phylo)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    #coord_fixed(ratio = max(HostMatrixdf$Space)/max(HostMatrixdf$Phylo))
    labs(title = "Space ~ Phylogeny"), #+ ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(title = "Sharing ~ Space"),# + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Phylogeny")# + ggpubr::stat_cor(method = "spearman")
  
), ncol = 3)


arrange_ggplot2(list(
  
  ggplot(HostMatrixdf[-LowerHosts,], aes(Space, PropVirus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Space") + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-LowerHosts,], aes(Space, PropVirus2)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing2 ~ Space") + ggpubr::stat_cor(method = "spearman")
  
), ncol = 2)


arrange_ggplot2(list(
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, PropVirus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Space") + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(HostMatrixdf[-HostThemselves,], aes(Space, PropVirus2)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing2 ~ Space") + ggpubr::stat_cor(method = "spearman")
  
), ncol = 2)


library(GGally)

HostThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

png(filename = "Figures/Pairwise Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
ggpairs(HostMatrixdf[-HostThemselves,], 
                lower = list(continuous = wrap("smooth", method = "gam"))) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space and Phylogeny Correlate with Viral Sharing")
dev.off()

# For Viruses ####

VirusThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(VirusAdj[FVN,FVN], diag = T) &
          lower.tri(VirusAdj[FVN,FVN], diag  = T))

png(filename = "Figures/Pairwise Virus Similarity Matrix Correlations.jpg", units = "mm", width = 200, height = 200, res = 300)
ggpairs(VirusHostMatrixdf[-VirusThemselves,], 
                lower = list(continuous = wrap("smooth", method = "gam"))) + # will take a while
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggtitle("Space Correlates with Host Sharing")
dev.off()

