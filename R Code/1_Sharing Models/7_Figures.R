
# Generating Figures ####

library(ggplot2); library(ggregplot); library(colorspace)

# 1.	Panel: Relationships/model outputs from: 
# Space and phylogeny; phylogeny and viral sharing; space and viral sharing for all data

load("Output Files/FitList.Rdata")

ggplot(FinalHostMatrix, aes(Space, Phylo)) + geom_point()

ggplot(FitList[[1]], aes(Space, Fit)) + 
  facet_wrap(~Phylo2) +
  geom_line()

FitList[[1]] %>% #filter(Space == mean(DataList[[1]]$Space)) %>% 
                        #DietSim == mean(DataList[[Resps[r]]]$DietSim)) %>% 
  ggplot(aes(Phylo2, Fit)) + geom_point()


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

load("Output Files/GridDegree2.Rdata")
load("Output Files/GridDegree4.Rdata")
load("Output Files/GridDegreeSum2.Rdata")
load("Output Files/GridDegreeSum4.Rdata")

PlotGrids <- GridDegree2

PlotGrids %>% filter(Metric == "AllPredDegree") %>% mutate(Degree = ifelse(Degree>320, 320, Degree)) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(AllPredDegree = "All Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  ggsave("Figures/All Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

PlotGrids %>% filter(Metric == "InDegree")  %>% mutate(Degree = ifelse(Degree>200, 200, ifelse(Degree<30,40,Degree))) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(InDegree = "Within-Order Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  ggsave("Figures/In Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

PlotGrids %>% filter(Metric == "OutDegree")  %>% mutate(Degree = ifelse(Degree>170, 170, ifelse(Degree<110,110,Degree))) %>%
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

# No. hosts versus predictability ####

ggplot(GAMValidSummary, aes(log10(NHosts), log10(MeanRank))) + geom_smooth() + geom_text(aes(label = Virus))

# Correlations among degree measures ####

GGally::ggpairs(Hosts %>% select(contains("Degree")), lower = list(continuous = "smooth"))

# Numbers of species versus centrality ####

Panth1 %>% group_by(hOrder) %>%
  summarise(Number = n(),
            AllPredDegree = mean(AllPredDegree),
            InDegree = mean(InDegree),
            OutDegree = mean(OutDegree)) %>%
  gather(key = "Metric", value = "Degree", contains("Degree"))

# Correlation between mean rank of focal host predictions and the proportion of links they're present for

ggplot(GAMValidSummary, aes(log10(MeanRank), MeanCount1)) + geom_point() + 
  geom_smooth(method = lm)

