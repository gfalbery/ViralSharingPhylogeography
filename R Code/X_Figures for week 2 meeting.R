
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


