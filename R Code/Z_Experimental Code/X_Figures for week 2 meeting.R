
# Figures for meeting with KO, NR, EE

library(grid)

# Are domestic viruses more central? ####

BarGraph(Hosts, "hDom", "Eigenvector", text = "N") +
  scale_fill_brewer(palette = AlberPalettes[4]) +
  theme(legend.position = "none") +
  labs(x = "Domestic", y = "Eigenvector", 
       title = "Domestic Animals are More Central") +
  ggsave("Figures/Domestic Centrality.jpeg", 
         units = "mm", height = 100, width = 100, dpi = 300)

# Yes, but sampling biases:

ggplot(Hosts, aes(log(hAllZACites+1), kader:::cuberoot(Eigenvector), colour = hDom)) + geom_point() + geom_smooth() +
  scale_fill_brewer(palette = AlberPalettes[4]) +
  labs(colour = "Domestic", title = "Publication Bias in Centrality") +
  ggsave("Figures/Domestic Publication Bias.jpeg", 
         units = "mm", height = 100, width = 100, dpi = 300)

# Space and Phylogeny correlate with viral sharing

jpeg("Figures/Space and phylogeny correlate with viral sharing.jpeg", units = "mm", width = 200, height = 100, res = 300)

arrange_ggplot2(list(
  ggplot(LongMatrixdf[-HostThemselves,], aes(Space, Phylo)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    #coord_fixed(ratio = max(LongMatrixdf$Space)/max(LongMatrixdf$Phylo))
    labs(title = "Space ~ Phylogeny"), #+ ggpubr::stat_cor(method = "spearman"),
  
  ggplot() + theme_void(),
  
  ggplot(LongMatrixdf[-HostThemselves,], aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2) +
    labs(title = "Sharing ~ Space"),# + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(LongMatrixdf[-HostThemselves,], aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    labs(title = "Sharing ~ Phylogeny")# + ggpubr::stat_cor(method = "spearman")
  
), ncol = 2)

dev.off()

cor1 <- with(LongMatrixdf[-HostThemselves,], cor.test(Space, Phylo, method = "spearman")$estimate) %>% round(2)
cor2 <- with(LongMatrixdf[-HostThemselves,], cor.test(Space, Virus, method = "spearman")$estimate) %>% round(2)
cor3 <- with(LongMatrixdf[-HostThemselves,], cor.test(Phylo, Virus, method = "spearman")$estimate) %>% round(2)

p1 <- with(LongMatrixdf[-HostThemselves,], cor.test(Space, Phylo, method = "spearman")$p.value)
p2 <- with(LongMatrixdf[-HostThemselves,], cor.test(Space, Virus, method = "spearman")$p.value)
p3 <- with(LongMatrixdf[-HostThemselves,], cor.test(Phylo, Virus, method = "spearman")$p.value)

ThemeRevert()

jpeg("Figures/Space and phylogeny correlate with viral sharing.jpeg", units = "mm", width = 150, height = 150, res = 300)

arrange_ggplot2(list(
  ggplot(LongMatrixdf[-HostThemselves,], aes(Space, Phylo)) + 
    geom_point(colour = AlberColours[1], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    #labs(title = "Phylogeny ~ Space") +
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") + 
    AlberTheme,
  
  ggplot(LongMatrixdf[1,], aes(Space, Space)) + 
    lims(x = c(0,1), y = c(0,1)) + theme_void() + 
    
    geom_text(data = data.frame(x = 0.5, y = 0.65), 
              inherit.aes = F, aes(x, y),
              label = paste0("R = ",cor1,", p < E-16"),
              colour = AlberColours[1],
              size = 5) + 
    
    geom_text(data = data.frame(x = 0.5, y = 0.5), 
              inherit.aes = F, aes(x, y),
              label = paste0("R = ",cor2,", p < E-16"),
              colour = AlberColours[2],
              size = 5) +
    
    geom_text(data = data.frame(x = 0.5, y = 0.35), 
              inherit.aes = F, aes(x, y),
              label = paste0("R = ",cor3,", p < E-16"),
              colour = AlberColours[3],
              size = 5),
  
  ggplot(LongMatrixdf[-HostThemselves,], aes(Space, Virus)) + 
    geom_point(colour = AlberColours[2], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    #labs(title = "Sharing ~ Space")  + 
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    AlberTheme,# + ggpubr::stat_cor(method = "spearman"),
  
  ggplot(LongMatrixdf[-HostThemselves,], aes(Phylo, Virus)) + 
    geom_point(colour = AlberColours[3], alpha = 0.3) + geom_smooth(colour = "black", fill = NA) +
    #labs(title = "Sharing ~ Phylogeny") + 
    stat_smooth(geom = "ribbon", fill = NA, lty = 2, col = "black") +
    AlberTheme# + ggpubr::stat_cor(method = "spearman")
  
), ncol = 2)

dev.off()

# Do domestics harbour more zoonoses?

jpeg("Figures/Domestics do not harbour more zoonoses.jpeg", units = "mm", width = 200, height = 100, res = 300)

arrange_ggplot2(list(
  
  BarGraph(Hosts, "Domestic", "hZoonosisProp", text = "N") + scale_fill_brewer(palette = AlberPalettes[4]) +
    theme(legend.position = "none"),
  
  BarGraph(Hosts, "Domestic", "hZoonosisCount", text = "N") + scale_fill_brewer(palette = AlberPalettes[4]) +
    theme(legend.position = "none")
  
), ncol = 2)

dev.off()

jpeg("Figures/Domestics Viruses are less likely to be zoonotic.jpeg", units = "mm", width = 200, height = 100, res = 300)

arrange_ggplot2(list(
  
  BarGraph(Viruses, "Domestic", "Human", text = "N") + scale_fill_brewer(palette = AlberPalettes[2]) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = 0:1),
  
  BarGraph(Viruses, "Domestic", "IsZoonotic", text = "N") + scale_fill_brewer(palette = AlberPalettes[2]) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = 0:1),
  
  BarGraph(Viruses, "Domestic", "IsZoonotic.stringent", text = "N") + scale_fill_brewer(palette = AlberPalettes[2]) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = 0:1)
  
), ncol = 3)

dev.off()

Viruses %>% BarGraph("Domestic", "Human", "Wildlife", text = "N") + scale_fill_brewer(palette = AlberPalettes[2]) +
  scale_x_continuous(breaks = 0:1) + ggsave("Figures/Domestic Virus_Wildlife Virus interactions.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

BarBarGraph(Viruses, "DomWild", "HostRangeMean", "Human") + 
  scale_color_manual(values = brewer.pal(6, AlberPalettes[4])[c(6,5,1)]) +
  ggsave("Figures/Domestic Virus_Wildlife Virus interactions2.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

BarBarGraph(Viruses, "DomWild", "Centroid_Human_Distance", "Human") + 
  scale_color_manual(values = brewer.pal(6, AlberPalettes[4])[c(6,5,1)])
ggsave("Figures/Domestic Virus_Wildlife Virus interactions3.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

jpeg("Figures/Domestics do not host more zoonoses2.jpeg", units = "mm", width = 200, height = 100, res = 300)

arrange_ggplot2(list(
  
  BarGraph(Hosts, "WellsDoms", "hZoonosisProp", text = "N") + scale_fill_brewer(palette = AlberPalettes[4]) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = 0:1),
  
  BarGraph(Hosts, "WellsDoms", "hZoonosisCount", text = "N") + scale_fill_brewer(palette = AlberPalettes[4]) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = 0:1)
  
), ncol = 2)

dev.off()

# Domestic Biases ####

ggplot(Hosts, aes(log(DOMYearBP), c(hZoonosisCount))) + geom_point() +
  labs(y = "Zoonosis Count") +
  ggrepel::geom_text_repel(aes(label = Sp), 
                           data = Hosts[Hosts$hZoonosisCount>10,],
                           arrow = arrow(length = unit(0.03, "npc")),
                           segment.size = 0.2, segment.color = "black") +
  ggsave("Figures/Time since domestication correlates with zoonosis count.jpeg",
         units = "mm", width = 100, height = 100, dpi = 300)

ggplot(Hosts, aes(log(DOMYearBP), c(hZoonosisProp))) + geom_point() +
  labs(y = "Proportion Zoonotic Viruses") +
  ggrepel::geom_text_repel(aes(label = Sp), 
                           data = Hosts[Hosts$hZoonosisProp>0.5,],
                           arrow = arrow(length = unit(0.03, "npc")),
                           segment.size = 0.2, segment.color = "black") +
  coord_fixed(ratio = 7) +
  ggsave("Figures/Time since domestication and proportion zoonoses.jpeg",
         units = "mm", width = 100, height = 100, dpi = 300)

ggplot(Hosts, aes(c(hZoonosisCount),c(hZoonosisProp))) + geom_point() + geom_smooth()
ggplot(Hosts, aes(c(hZoonosisCount),c(hZoonosisProp))) + geom_point() + geom_smooth(method = lm)

ggplot(Hosts, aes(log(hAllZACites),c(hZoonosisProp))) + geom_point() + 
  geom_smooth(method = lm, formula = y~poly(x,2))

jpeg("Proportional spatial overlap with other mammals explains zoonosis count.Rdata")
ggpairs(Hosts[,c(paste0("S",c("",c(2,4,5,8,10)*10)),"S.Greg1","S.Greg2","hZoonosisCount", "hZoonosisProp")] %>% 
          mutate(lZoos = log(hZoonosisCount + 1)), 
        lower = list(continuous = wrap("smooth", method = "loess"))) + 
  AlberTheme
dev.off()