
# Figures for Preliminary Model Output ####

#devtools::install_github("gfalbery/ggregplot")
library(ggregplot); library(ggplot2); library(RColorBrewer)

ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

AlberTheme <- theme_bw() +
  theme(axis.title.x = element_text(vjust = -0.35, 
                                    size = 12, 
                                    colour = "black"), 
        axis.title.y = element_text(vjust = 1.2, 
                                    size = 12, 
                                    colour = "black"),
        strip.background = element_rect(fill = "white", colour = "dark grey"))

theme_set(AlberTheme)


# ZI Model Output ####

ZISols <- summary(ZIModel)$solutions %>% as.data.frame %>% mutate(
  Component = rep(c("Count", "ZI"), dim(ZISols)[1]/2),
  Variable = rep(c("Intercept", "Space", "Phylogeny", "Citations", "DomDom", "DomWild", 
                   "Phylo:Space"),
                 each = 2),
  Name = paste(Component, Variable, sep = ":")
) %>% rename(Lower = "l-95% CI", Upper = "u-95% CI", Estimate = "post.mean")

ZISols$Estimate[ZISols$Component=="ZI"] <- -ZISols$Estimate[ZISols$Component=="ZI"]
ZISols[ZISols$Component=="ZI", c("Lower", "Upper")] <- -ZISols[ZISols$Component=="ZI", c("Upper","Lower")]

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = Component)) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + coord_flip() + AlberTheme +
  ggtitle("Zero-Inflated Model Output (no spatial zeroes)") +
  ggsave("Figures/Zero-Inflated Model Output (no spatial zeroes).jpeg", units = "mm", width = 150, height = 150, dpi = 300)

# Space and phylogeny correlate with viruses ####

jpeg("Figures/Space and Phylogeny correlate with viral sharing with no spatial zeroes.jpeg",
     units = "mm", height = 150, width = 150, res = 300)

list(
  
  HostMatrixdf %>% filter(!Sp==Sp2) %>%  ggplot(aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, colour = "black", method = lm) +
    labs(x = "Space Shared", y = "Viruses Shared", title = "Virus ~ Space"),
  
  HostMatrixdf %>% filter(Space>0, !Sp==Sp2) %>%  ggplot(aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, colour = "black", method = lm) +
    lims(x = c(0,1)) +
    labs(x = "Space Shared", y = "Viruses Shared", title = "Virus ~ Space, >0 Space"),
  
  HostMatrixdf %>% filter(!Sp==Sp2) %>%  ggplot(aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, colour = "black", method = lm) +
    labs(x = "Genetic Similarity", y = "Viruses Shared",title = "Virus ~ Phylogeny"),
  
  HostMatrixdf %>% filter(Space>0, !Sp==Sp2) %>%  ggplot(aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, colour = "black", method = lm) +
    lims(x = c(0,1)) +
    labs(x = "Genetic Similarity", y = "Viruses Shared", title = "Virus ~ Phylogeny, >0 Space")
  
) %>% arrange_ggplot2(nrow = 2)

dev.off()

# Looking at different slopes of phylo:virus ####

HostMatrixdf %>% filter(!Sp == Sp2) %>% ggplot(aes(Phylo, Virus, colour = SpaceQuantile)) + 
  geom_point(alpha = 0.3) + geom_smooth(method = lm) + facet_wrap(~SpaceQuantile) +
  ggtitle("Space sharing changes the slope of Virus ~ Phylogeny") +
  ggsave("Figures/Phylogeny:Virus slope varies with space.jpeg", units = "mm", height = 150, width = 200, dpi = 300)

# Second part: centrality analyses ####

jpeg("Figures/Different components of centrality correspond differently to space and phylogeny.jpeg", units = "mm", width = 150, height = 10, res = 300)

lapply(CentralityList, function(a) INLADICFig(a[1:4], ModelNames = ModelNames) + labs(x = "Model Spec")) %>% 
  arrange_ggplot2(nrow = 3)

dev.off()

jpeg("Figures/Space and phylogeny do not strongly affect centrality effect estimates.jpeg", units = "mm", width = 300, height = 150, res = 300)

lapply(CentralityList, function(a) Efxplot(a[1:4], ModelNames = ModelNames) + theme(legend.position = "top") + labs(colour = NULL)) %>% 
  arrange_ggplot2(ncol = 3)

dev.off()

ggField(CentralityList[[3]][[2]], WorldMesh) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long/50000, lat/50000, group = group)) +
  geom_point(data = TestHosts[,c("LongMean", "LatMean")], aes(LongMean, LatMean), inherit.aes = F) + 
  scale_fill_brewer(palette = AlberPalettes[2]) + 
  ggtitle("Spatial Distribution of Eigenvector Centrality") +
  ggsave("Figures/Eigenvector Centrality Spatial Distribution.jpeg", units = "mm", height = 150, width = 200, dpi = 300)


# Model outputs ####

jpeg("Figures/DIC from centrality models excluding hOrder.jpeg", units = "mm", width = 250, height = 200, res = 300)
lapply(CentralityList, function(a) INLADICFig(a, ModelNames = c(ModelNames,"SPDE+SMat")) + 
         theme(legend.position = "none") + labs(x  = "Model")
       ) %>% #)) %>% 
  arrange_ggplot2(nrow = 3)
dev.off()

jpeg("Figures/DIC from centrality models including hOrder.jpeg", units = "mm", width = 250, height = 200, res = 300)
lapply(SaveCentralityList, function(a) INLADICFig(a, ModelNames = c(ModelNames,"SPDE+SMat")) + 
         theme(legend.position = "none") + labs(x  = "Model")
       ) %>% #, ModelNames = ModelNames)) %>% 
  arrange_ggplot2(nrow = 3)
dev.off()

for(i in 1:length(SaveCentralityList)){
  Efxplot(SaveCentralityList[[i]], ModelNames = c(ModelNames,"SPDE+SMat")) + 
    ggtitle(paste0(names(SaveCentralityList)[i], " Model Output")) +
    ggsave(paste0("Figures/Model Outputs including hOrder_",names(SaveCentralityList)[i],".jpeg"), units = "mm", width = 200, height = 200, dpi = 300)
}



