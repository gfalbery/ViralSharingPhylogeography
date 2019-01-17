# Viral Distances to Humans and Domestic Animals ####

library(igraph); library(ggregplot)

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

vHumanDist <- distances(bipgraph, v = V(bipgraph)[1:(dim(Viruses)[1])],
                        to = V(bipgraph)[which(names(V(bipgraph)) == "Homo_sapiens")])

vHumanDist[,1] <- ifelse(vHumanDist[,1] == Inf, Inf, (vHumanDist[,1]-1)/2)

vHumanDist2 <- data.frame(Sp = rownames(vHumanDist),
                          HumanDist = vHumanDist[,1],
                          Domestic = Viruses$Domestic,
                          Family = Viruses$vFamily)

vDomDist <- distances(bipgraph, v = V(bipgraph)[1:(dim(Viruses)[1])],
                      to = V(bipgraph)[which(names(V(bipgraph)) %in% Domestics)])

vDomDist <- (vDomDist-1)/2

vDomDist2 <- data.frame(Sp = rownames(vDomDist),
                        MeanDomDist = apply(vDomDist, 1, function(a) mean(a[!a=="Inf"], na.rm = T)),
                        MinDomDist = apply(vDomDist, 1, function(a) ifelse(all(a==Inf),Inf,min(a[!a==Inf], na.rm = T))),
                        MaxDomDist = apply(vDomDist, 1, function(a) max(a, na.rm = T)),
                        Human = Viruses$Human,
                        Family = Viruses$vFamily)

# Merging distances and traits ####

Viruses <- merge(Viruses, 
                 vDomDist2[,c("Sp", "MeanDomDist", "MinDomDist", "MaxDomDist")], 
                 by = "Sp", all.x = T)

Viruses <- merge(Viruses, 
                 vHumanDist2[,c("Sp", "HumanDist")], 
                 by = "Sp", all.x = T)

# Making long dataset with all links ####

vDomDistLong <- reshape2::melt(vDomDist)

ggplot(vDomDistLong, aes(Var1, value, colour = Var1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  geom_point(position = position_jitter(w = 0.5)) 

# Attempting some plots ####

# Distance to Domestics ####

BarGraph(vDomDist2[!NARows(vDomDist2),], "Human", "Min", text = "N") + 
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = c(0:1)) +
  coord_fixed(ratio = 1/2) +
  labs(x =  "Infects Humans", 
       y = "Min Distance to Domestics",
       title = "vHuman Min Distance to Domestics")  +
  ggsave("Figures/vHuman Min Distance to Domestics.tiff", units = "mm", width = 100, height = 100, dpi = 300)

BarGraph(vDomDist2[!NARows(vDomDist2),], "Human", "Mean", text = "N") + 
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = c(0:1)) +
  coord_fixed(ratio = 1/2) +
  labs(x =  "Infects Humans", 
       y = "Mean Distance to Domestics",
       title = "vHuman Mean Distance to Domestics")  +
  ggsave("Figures/vHuman Mean Distance to Domestics.tiff", units = "mm", width = 100, height = 100, dpi = 300)

BarGraph(vDomDist2[!NARows(vDomDist2),], "Family", "Min", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_fixed(ratio = nunique(vDomDist2$Family)) +
  labs(x =  "vFamily", 
       y = "Min Distance from Domestics",
       title = "vFamily Min Distance to Domestics")  +
  ggsave("Figures/vFamily Min Distance to Domestics.tiff", units = "mm", width = 200, height = 200, dpi = 300)

BarGraph(vDomDist2[!NARows(vDomDist2),], "Family", "Mean", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_fixed(ratio = nunique(vDomDist2$Family)/5) +
  labs(x =  "vFamily", 
       y = "Mean Distance from Domestics",
       title = "vFamily Mean Distance to Domestics")  +
  ggsave("Figures/vFamily Mean Distance to Domestics.tiff", units = "mm", width = 200, height = 200, dpi = 300)

# Distance from Zoonosis ####

BarGraph(vHumanDist2[!NARows(vHumanDist2),], "Domestic", "HumanDist", text = "N") + 
  theme(legend.position = "none") + 
  labs(x =  "Infects Domestics", 
       y = "Distance to Humans",
       title = "vDomestic Distance to Humans") +
  scale_x_continuous(breaks = c(0:1)) + 
  coord_fixed(ratio = 2/3) +
  ggsave("Figures/vDomestic Distance to Humans.tiff", units = "mm", width = 150, height = 150, dpi = 300)

BarGraph(vHumanDist2[!NARows(vHumanDist2),], "Family", "HumanDist", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Family", 
       y = "Min Distance to Humans",
       title = "vFamily Distance to Humans") + 
  coord_fixed(ratio = nunique(Viruses$vFamily)/3) +
  theme(legend.position = "none") +
  ggsave("Figures/vFamily Distance to Humans.tiff", units = "mm", width = 150, height = 150, dpi = 300)

ggMMplot(vHumanDist2, "Domestic", "HumanDist") +
  scale_fill_discrete(limits = c(0:2,"Inf")) + 
  labs(x = "Infects Domestics", 
       y = "vDomestic Distance from Humans",
       title = "Domestic Distances to Humans") + coord_fixed() +
  ggsave("Figures/vDomestic Distance to Humans Mosaic.tiff", units = "mm", width = 100, height = 100, dpi = 300)

ggMMplot(vHumanDist2, "Family", "HumanDist") +
  scale_fill_discrete(limits = c(0:2,"Inf")) + 
  labs(x = "Infects Domestics", 
       y = "vFamily Distance from Humans",
       title = "Family Distances to Humans") + coord_fixed() +
  ggsave("Figures/vFamily Distance to Humans Mosaic.tiff", units = "mm", width = 100, height = 100, dpi = 300)

