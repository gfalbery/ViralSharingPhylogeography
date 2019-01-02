# Distances ####

NARows <-function(df, vars){
  apply(df[,vars], 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

vHumanDist <- distances(bipgraph, v = V(bipgraph)[1:(dim(Viruses)[1])],
                        to = V(bipgraph)[which(names(V(bipgraph)) == "Homo_sapiens")], weights=NA)

vDomDist <- distances(bipgraph, v = V(bipgraph)[1:(dim(Viruses)[1])],
                        to = V(bipgraph)[which(names(V(bipgraph)) %in% Domestics)], weights=NA)

vDomDist2 <- data.frame(Sp = rownames(vDomDist),
                        Mean = apply(vDomDist, 1, function(a) mean(a[!a=="Inf"], na.rm = T)),
                        Min = apply(vDomDist, 1, function(a) min(a[!a=="Inf"], na.rm = T)),
                        Max = apply(vDomDist, 1, function(a) max(a, na.rm = T)),
                        Human = Viruses$Human,
                        Family = Viruses$vFamily)

BarGraph(vDomDist2, "Human", "Min", text = "N") + 
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = c(0:1)) +
  coord_fixed(ratio = 2) +
  labs(x =  "Infects Humans", 
       y = "Min Distance to Domestics",
       title = "vFamily Min Distance to Domestics")  +
  ggsave("Figures/vHuman Min Distance to Domestics.tiff", units = "mm", width = 100, height = 100, dpi = 300)

BarGraph(vDomDist2, "Family", "Min", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  coord_fixed(ratio = nunique(vDomDist2$Family)) +
  
  
  labs(x =  "vFamily", 
       y = "Min Distance from Domestics",
       title = "vFamily Min Distance to Domestics")  +
  ggsave("Figures/vFamily Min Distance to Domestics.tiff", units = "mm", width = 200, height = 200, dpi = 300)


vHumanDist2 <- vHumanDist[!NARows(vHumanDist),]

vHumanDist2 <- data.frame(Sp = rownames(vHumanDist2),
                          Mean = apply(vHumanDist2, 1, function(a) mean(a, na.rm = T)),
                          Min = apply(vHumanDist2, 1, function(a) min(a, na.rm = T)),
                          Max = apply(vHumanDist2, 1, function(a) max(a, na.rm = T)),
                          Domestic = Viruses$Domestic[!NARows(vHumanDist)],
                          Family = Viruses$vFamily[!NARows(vHumanDist)])

BarGraph(vHumanDist2, "Domestic", "Mean") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x =  "Infects Domestics", 
       y = "Mean Distance to Humans") +
  ggsave("Figures/vDomestic Mean Distance to Humans.tiff", units = "mm", width = 150, height = 150, dpi = 300)

BarGraph(vHumanDist2, "Family", "Min", text = "N") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Family", 
       y = "Mean Distance to Humans",
       title = "vFamily Min Distance to Humans") + 
  coord_fixed(ratio = nunique(Viruses$vFamily)/max(vHumanDist2$Mean)) +
  theme(legend.position = "none") +
  ggsave("Figures/vFamily Min Distance to Humans.tiff", units = "mm", width = 150, height = 150, dpi = 300)

ggMMplot(vHumanDist2, "Domestic", "Min") +
  scale_fill_discrete(limits = 0:3) + 
  labs(x = "Infects Domestics", 
       y = "Min Distance from Humans",
       title = "Domestic Min Distances to Humans") + coord_fixed() +
  ggsave("Figures/vDomestic Min Distance to Humans.tiff", units = "mm", width = 100, height = 100, dpi = 300)

ggMMplot(vHumanDist2, "Domestic", "Max") +
  scale_fill_discrete(limits = 0:3) + 
  labs(x = "Infects Domestics", 
       y = "Max Distance from Humans",
       title = "Domestic Max Distances to Humans") + coord_fixed() +
  ggsave("Figures/vDomestic Max Distance to Humans.tiff", units = "mm", width = 100, height = 100, dpi = 300)

# Putting in long format ####

vHumanDistLong <- reshape2::melt(vHumanDist)
vHumanDistLong$Domestic <- Viruses$Domestic

ggplot(vHumanDistLong, aes(as.factor(Domestic), value, colour = as.factor(Domestic))) + 
  ggforce::geom_sina() + 
  labs(x = "Infects Domestics", 
       y = "Mean Distances from Humans",
       title = "Domestic Distances to Humans") + coord_fixed(ratio = 2/3) +
  theme(legend.position = "none") +
  ggsave("Figures/Domestic-Human Total Distances 1.tiff", units = "mm", width = 100, height = 100, dpi = 300)

ggMMplot(vHumanDistLong, "Domestic", "value") +
  ggtitle("Domestic Distance to Humans") + coord_fixed() +
  ggsave("Figures/Domestic-Human Total Distances 2.tiff", units = "mm", width = 100, height = 100, dpi = 300)


# Distances from other types of Hosts ####
# dist.from.NYT # in case trying to google/evernote search the original code

hDomDist <- distances(Hostgraph, v = V(Hostgraph),
                      to = V(Hostgraph)[Hosts$hDom == "domestic"], weights=NA)

hHumanDist <- distances(Hostgraph, v = V(Hostgraph),
                        to = V(Hostgraph)[Hosts$hDom == "human"], weights=NA)

hHumanDist2 <- as.data.frame(hHumanDist)
hHumanDist2$hDom <- Hosts$hDom

hHumanDist2 <- hHumanDist2[!hHumanDist2$Homo_sapiens %in% c("Inf",0),]
BarGraph(hHumanDist2, "hDom", "Homo_sapiens")



