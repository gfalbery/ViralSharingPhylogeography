
# Generating Figures ####

library(ggplot2); library(ggregplot); library(colorspace)

# 1.	Panel: Relationships/model outputs from: 
# Space and phylogeny; phylogeny and viral sharing; space and viral sharing for all data

ggplot(FinalHostMatrix, aes(Space, Phylo)) + geom_point()



# 2.	Predicted degree centrality with and without random effect for ~
# 700 known viral hosts (possibly supplemental)

ggplot(Hosts, aes(Degree, PredDegree1)) + geom_point() +
  geom_smooth()

ggplot(Hosts, aes(Degree, PredDegree1b)) + geom_point() +
  geom_smooth()

ggplot(Hosts, aes(Degree, AllPredDegree)) + geom_point() +
  geom_smooth()

# 3.	Mammal order level centrality (bar plots)


BarGraph(Panth1, "hOrder", "AllPredDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  labs(x = "Order", y = "Degree Centrality", title = "All Links") +  
  scale_fill_discrete_sequential(palette = AlberPalettes[1]) +
  ggsave("Figures/AllPredDegree.jpeg", units = "mm", height= 120, width = 200)

BarGraph(Panth1, "hOrder", "InDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  labs(x = "Order", y = "Degree Centrality", title = "Within-Order Links") +  
  scale_fill_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("Figures/InDegree.jpeg", units = "mm", height= 120, width = 200)

BarGraph(Panth1, "hOrder", "OutDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  labs(x = "Order", y = "Degree Centrality", title = "Out-of-Order Links") +  
  scale_fill_discrete_sequential(palette = AlberPalettes[3]) +
  ggsave("Figures/OutDegree.jpeg", units = "mm", height= 120, width = 200)




# 4.	Gridcell level centrality measure (map), to see if overlaps w species diversity (optional)

load("Output Files/GridDegree.Rdata")

GridDegree2 %>% filter(Metric == "AllPredDegree") %>% filter(Degree<310) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(AllPredDegree = "All Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  ggsave("Figures/All Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

GridDegree2 %>% filter(Metric == "InDegree") %>% filter(Degree<180) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(InDegree = "Within-Order Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  ggsave("Figures/In Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

GridDegree2 %>% filter(Metric == "OutDegree") %>% filter(Degree<190) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(OutDegree = "Out-of-Order Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[3]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[3]) +
  ggsave("Figures/Out Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)



# 5.	Panel: viral trait plots (RNA vs. DNA; vector-borne or not; etc)


# 6.	Panel of maps: Illustrative maps of host ranges for mammal species from host predictions (currently unidentified) for specific viruses, e.g. CCHFV; Nipah virus; Ebola; etc. (supplement will list out probability of being a host for all viruses)

PredHostPlot("Crimean-congo_haemorrhagic_fever_virus", focal = c(1,0), facet = T)

list(PredHostPlot("Andes_virus", focal = 0), PredHostPlot("Andes_virus", focal = 1)) %>% 
  arrange_ggplot2




# Supplementary? #####

# No. hosts versus predictability
ggplot(ValidSummary, aes(log10(NHosts), log10(MeanRank))) + geom_smooth() + geom_text(aes(label = Virus))
ggplot(GAMValidSummary, aes(log10(NHosts), log10(MeanRank))) + geom_smooth() + geom_text(aes(label = Virus))






