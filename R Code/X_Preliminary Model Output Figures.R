
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

Efxplot(list(ZIModelNoSpace)) + 
  labs()


# Space and phylogeny correlate with viruses ####

jpeg("Figures/Space and Phylogeny correlate with viral sharing with no spatial zeroes.jpeg",
     units = "mm", height = 150, width = 150, res = 300)

list(
  
  HostMatrixdf %>% filter(!Sp==Sp2) %>%  ggplot(aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(x = "Space Shared", y = "Viruses Shared", title = "Virus ~ Space"),
  
  HostMatrixdf %>% filter(Space>0, !Sp==Sp2) %>%  ggplot(aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    lims(x = c(0,1)) +
    labs(x = "Space Shared", y = "Viruses Shared", title = "Virus ~ Space, >0 Space"),
  
  HostMatrixdf %>% filter(!Sp==Sp2) %>%  ggplot(aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(x = "Genetic Similarity", y = "Viruses Shared",title = "Virus ~ Phylogeny"),
  
  HostMatrixdf %>% filter(Space>0, !Sp==Sp2) %>%  ggplot(aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(method = lm, colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
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








