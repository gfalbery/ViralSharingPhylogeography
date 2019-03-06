
# Generating Figures ####

library(ggplot2); library(ggregplot); library(colorspace)

# 1.	Panel: Relationships/model outputs from: 
# Space and phylogeny; phylogeny and viral sharing; space and viral sharing for all data

load("Output Files/FitList.Rdata")

ggplot(FinalHostMatrix, aes(Space, Phylo)) + 
  geom_point(alpha = 0.1, colour = AlberColours[1]) +
  coord_fixed() +
  geom_smooth(colour = "black", se = F) +
  ggsave("Figures/Space_Phylo.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

FitList[["VirusBinary"]] %>% 
  filter(Phylo == last(PhyloRange)) %>%
  filter(DietSim == last(DietRange)) %>%
  ggplot(aes(Space, Fit)) + geom_line() +
  lims(y = c(0,1)) +
  labs(y = "Predicted Viral Sharing", x = "Geographic Overlap") +
  ggsave("Figures/Phylo Predictions.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

FitList[["VirusBinary"]] %>% 
  filter(Space == last(SpaceRange)) %>%
  filter(DietSim == last(DietRange)) %>%
  ggplot(aes(Phylo, Fit)) + geom_line() +
  lims(y = c(0,1)) +
  labs(y = "Predicted Viral Sharing", x = "Phylogenetic Similarity") +
  ggsave("Figures/Space Predictions.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

jpeg("Figures/Model Predictions.jpeg", units = "mm", width = 300, height = 150, res = 300)

list(
  FitList[["VirusBinary"]] %>% filter(!Space == last(unique(FitList[["VirusBinary"]]$Space))&
                                        DietSim == last(unique(FitList[["VirusBinary"]]$DietSim))) %>%
    
    ggplot(aes(Phylo, Fit, colour = Space)) + 
    geom_line(aes(group = as.factor(Space)), alpha = 0.3) +
    labs(y = "Predicted Viral Sharing", x = "Phylogenetic Similarity") +
    geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
  
  FitList[["VirusBinary"]] %>% filter(!Phylo == last(unique(FitList[["VirusBinary"]]$Phylo))&
                                        DietSim == last(unique(FitList[["VirusBinary"]]$DietSim))) %>%
    
    ggplot(aes(Space, Fit, colour = Phylo)) + 
    geom_line(aes(group = as.factor(Phylo)), alpha = 0.3) +
    labs(y = "Predicted Viral Sharing", x = "Geographic Overlap") +
    geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01)
  
) %>% arrange_ggplot2

dev.off()

# 2.	Predicted degree centrality with and without random effect for ~
# 700 known viral hosts (possibly supplemental)

ggplot(Hosts, aes(Degree, PredDegree1)) + geom_point() +
  geom_smooth()

ggplot(Hosts, aes(Degree, PredDegree1b)) + geom_point() +
  geom_smooth()

ggplot(Hosts, aes(Degree, AllPredDegree)) + geom_point() +
  geom_smooth()

# 3.	Mammal order level centrality (bar plots)

load("Output Files/Panth1.Rdata")

BarGraph(Panth1, "hOrder", "AllPredDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  labs(x = "Order", y = "Degree Centrality", title = "All Links") +  
  scale_fill_discrete_sequential(palette = AlberPalettes[1]) +
  ggsave("Figures/AllPredDegree.jpeg", units = "mm", height = 120, width = 200)

BarGraph(Panth1, "hOrder", "InDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  labs(x = "Order", y = "Degree Centrality", title = "Within-Order Links") +  
  scale_fill_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("Figures/InDegree.jpeg", units = "mm", height = 120, width = 200)

BarGraph(Panth1, "hOrder", "OutDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  labs(x = "Order", y = "Degree Centrality", title = "Out-of-Order Links") +  
  scale_fill_discrete_sequential(palette = AlberPalettes[3]) +
  ggsave("Figures/OutDegree.jpeg", units = "mm", height = 120, width = 200)




# 4.	Gridcell level centrality measure (map), to see if overlaps w species diversity (optional)

load("Output Files/GridDegree2.Rdata")
load("Output Files/GridDegree4.Rdata")
load("Output Files/GridDegreeSum2.Rdata")
load("Output Files/GridDegreeSum4.Rdata")

PlotGrids <- GridDegree2

PlotGrids %>% filter(Metric == "AllPredDegree") %>% mutate(Degree = ifelse(Degree>320, 320, Degree)) %>%
  ggplot(aes(x, y, fill = Degree)) + #, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(AllPredDegree = "All Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  ggsave("Figures/All Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

PlotGrids %>% filter(Metric == "InDegree")  %>% mutate(Degree = ifelse(Degree>200, 200, ifelse(Degree<30,40,Degree))) %>%
  ggplot(aes(x, y, fill = Degree)) + # , colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(InDegree = "Within-Order Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  ggsave("Figures/In Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

PlotGrids %>% filter(Metric == "OutDegree")  %>% mutate(Degree = ifelse(Degree>170, 170, ifelse(Degree<110,110,Degree))) %>%
  ggplot(aes(x, y, fill = Degree)) + #, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(OutDegree = "Out-of-Order Links"))) +
  coord_fixed() +  
  lims(x = c(80, 720)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[3]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[3]) +
  ggsave("Figures/Out Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

# 5.	Panel: viral trait plots (RNA vs. DNA; vector-borne or not; etc)



# 6.	Panel of maps: Illustrative maps of host ranges for mammal species from host predictions (currently unidentified) for specific viruses, e.g. CCHFV; Nipah virus; Ebola; etc. (supplement will list out probability of being a host for all viruses)

PredPlot(HostList = VirusAssocs[["Crimean-Congo_hemorrhagic_fever_virus"]], 
         Focal = c(1), 
         Facet = F,
         Validate = F) + 
  theme_void() + theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("CCHF Known Hosts") +
  ggsave("Figures/CCHF_Known.jpeg", units = "mm", width = 250, height = 150)

PredPlot(HostList = VirusAssocs[["Crimean-Congo_hemorrhagic_fever_virus"]], 
         Focal = c(0), 
         Facet = F,
         Validate = F, Threshold = 15) + 
  theme_void() + theme(plot.title = element_text(hjust=0.5)) +
  #ggtitle("CCHF Model Success") +
  ggsave("Figures/CCHF_Predicted.jpeg", units = "mm", width = 250, height = 150)

PredPlot(HostList = VirusAssocs[["Crimean-Congo_hemorrhagic_fever_virus"]], 
         Focal = c(1,0), 
         Facet = T,
         Validate = T)[[2]] + 
  ggtitle("CCHF Model Success") +
  ggsave("Figures/CCHF_Success.jpeg", units = "mm", width = 100, height = 100)

PredPlot(Virus = "Crimean-Congo_hemorrhagic_fever_virus", 
         Focal = c(0), 
         Validate = F,
         Facet = F)

list(PredHostPlot("Andes_virus", focal = 0), PredHostPlot("Andes_virus", focal = 1)) %>% 
  arrange_ggplot2




# Supplementary? #####

# Model description: posterior draws ####

ggplot(PostList[["VirusBinary"]]$Space, aes(i, Fit, colour = Draw)) + geom_line(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "Space", y = "Model Estimate", title = "Posterior Draw Estimates") +
  geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01) +
  scale_colour_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("SIFigures/GAMPosteriors_Space.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

ggplot(PostList[["VirusBinary"]]$Phylo, aes(i, Fit, colour = Draw)) + geom_line(alpha = 0.3) + theme(legend.position = "none") +
  labs(x = "Phylo", y = "Model Estimate", title = "Posterior Draw Estimates") +
  geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01) +
  scale_colour_discrete_sequential(palette = AlberPalettes[1]) +
  ggsave("SIFigures/GAMPosteriors_Phylo.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

# Model description: data draws ####

ggplot(DrawList[["VirusBinary"]]$Space, 
       aes(i, Fit, colour = Iteration)) + 
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01) +
  labs(x = "Space", y = "Model Estimate", title = "Data Draw Estimates") +
  scale_colour_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("SIFigures/GAMDataDraws_Space.jpeg", 
         units = "mm", width = 100, height = 100, dpi = 300)

ggplot(DrawList[["VirusBinary"]]$Phylo, 
       aes(i, Fit, colour = Iteration)) + 
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01) +
  labs(x = "Phylo", y = "Model Estimate", title = "Data Draw Estimates") +
  scale_colour_discrete_sequential(palette = AlberPalettes[1]) +
  ggsave("SIFigures/GAMDataDraws_Phylo.jpeg", 
         units = "mm", width = 100, height = 100, dpi = 300)

# Species-level Correlations among degree predictions ####

jpeg("SIFigures/DegreePairs.jpeg", units = "mm", width = 200, height = 200, res = 300)

GGally::ggpairs(Hosts %>% select(contains("Degree")), lower = list(continuous = "smooth"))

dev.off()

# Order-level Correlations among degree predictions ####

BarBarGraph(Hosts, "hOrder", "Degree", "AllPredDegree") +
  ggtitle("Order-Level Predictions Do Not Correlate") +
  ggsave("SIFigures/OrderDegreePredicts.jpeg", units = "mm", height = 120, width = 125, dpi = 300)

# Numbers of species versus centrality ####

Panth1 %>% group_by(hOrder) %>%
  summarise(Number = n(),
            AllPredDegree = mean(AllPredDegree),
            InDegree = mean(InDegree),
            OutDegree = mean(OutDegree)) %>%
  gather(key = "Metric", value = "Degree", contains("Degree")) %>%
  ggplot(aes(Number, Degree)) + 
  ggtitle("Order Size") + #coord_fixed() +
  geom_point() + geom_smooth(method = lm, colour = "black") + facet_wrap(~Metric) +
  labs(x = "Number of Mammals") +
  ggsave("SIFigures/NHosts_Degree_Order.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

Panth1 %>% group_by(hOrder) %>%
  summarise(Number = n(),
            AllPredDegree = mean(AllPredDegree),
            InDegree = mean(InDegree),
            OutDegree = mean(OutDegree)) %>%
  gather(key = "Metric", value = "Degree", contains("Degree")) %>%
  ggplot(aes(log(Number), log(Degree))) + 
  ggtitle("Order Size") + #coord_fixed() +
  geom_point() + geom_smooth(method = lm, colour = "black") + facet_wrap(~Metric) +
  labs(x = "Number of Mammals") +
  ggsave("SIFigures/log_NHosts_Degree_Order.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

Panth1 %>% group_by(hOrder) %>%
  summarise(Number = n(),
            AllPredDegree = mean(AllPredDegree),
            InDegree = mean(InDegree),
            OutDegree = mean(OutDegree)) %>% lm(log(InDegree+1) ~ log(Number+1), data = .)

OrderLevelLinks %>% ggplot(aes(log(HostNumber), log(Degree+1))) + geom_point() + 
  coord_fixed() + 
  geom_smooth(colour = "black", fill = NA) + 
  stat_smooth(fill = NA, geom = "ribbon", lty = 2, colour = "black") +
  facet_wrap(~Metric, labeller = labeller(Metric = c("AllPredDegree" = "All Links",
                                                     "OutDegree" = "Out-of-Order Links",
                                                     "InDegree" = "Within-Order Links")))  +
  ggsave("SIFigures/log_NHosts_Degree_Order.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

hComboList %>% group_by(Iteration, Group) %>%
  summarise(HostNo = n(),
            Degree = sum(Degree)) %>% spread(value = c("HostNo"), key = c("Group")) %>%
  group_by(Iteration) %>%
  summarise(`1` = sum(`1`, na.rm = T), 
            `2` = sum(`2`, na.rm = T),
            Degree = mean(Degree)) %>% bind_cols(Combos) %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1))) + 
  geom_point() +
  labs(x = "log(Order 1 Hosts) + log(Order 2 Hosts)", 
       y = "log(Predicted Links + 1)",
       title = "Scaling of between-order links") +
  ggsave("SIFigures/BetweenOrderScaling.jpeg", units = "mm", height = 100, width = 100, dpi = 300)

# No. hosts versus predictability ####

ggplot(ValidSummary, aes(log10(NHosts), log10(MeanRank))) + 
  geom_point() +
  geom_smooth(colour = "black", fill = NA) + 
  stat_smooth(fill = NA, geom = "ribbon", lty = 2, colour = "black") +
  #geom_text(aes(label = Virus)) +
  labs(x = "log10(Number of Hosts)", y = "log10(Mean Rank of Focal Host)") +
  ggsave("SIFigures/HostNo_Prediction.jpeg", units = "mm", height = 100, width = 100, dpi = 300)

# Correlation between mean rank of focal host predictions and the proportion of links they're present for

ggplot(ValidSummary, aes(log10(MeanRank), MeanCount1)) + geom_point() + 
  geom_smooth(method = lm) +
  ggsave("SIFigures/Rank_Count.jpeg", units = "mm", height = 100, width = 100)

# RNA Viruses are harder to predict ####

BarGraph(ValidSummary, "vDNAoRNA", "MeanRank", text = "N") +
  scale_fill_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  theme(legend.position = "none") +
  ggsave("SIFigures/RNA_DNA_Prediction.jpeg", units = "mm", height = 100, width = 100)

# No other viral traits matter ####

jpeg("SIFigures/ViralTraits_Predictability.jpeg", units = "mm", width = 200, height = 200, res = 300)

VirusCovar %>% 
  lapply(function(a) BarGraph(ValidSummary, a, "MeanRank", text = "N")) %>% 
  arrange_ggplot2(ncol = 3)

dev.off()



