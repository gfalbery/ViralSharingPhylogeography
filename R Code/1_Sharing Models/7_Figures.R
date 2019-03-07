
# Generating Figures ####

library(ggplot2); library(ggregplot); library(colorspace)

# Figure 1_ space and phylogeny ####

load("Output Files/FitList.Rdata")

ggplot(FinalHostMatrix, aes(Space, Phylo)) + 
  geom_point(alpha = 0.1, colour = AlberColours[1]) +
  coord_fixed() +
  geom_smooth(colour = "black", se = F) +
  ggsave("Figures/Space_Phylo.jpeg", units = "mm", width = 100, height = 100, dpi = 300)

jpeg("Figures/Model Predictions.jpeg", units = "mm", width = 300, height = 150, res = 600)

list(
  FitList[["VirusBinary"]] %>% filter(!Space == last(unique(FitList[["VirusBinary"]]$Space))&
                                        DietSim == last(unique(FitList[["VirusBinary"]]$DietSim))) %>%
    
    ggplot(aes(Phylo, Fit, colour = Space)) + 
    geom_line(aes(group = as.factor(Space)), alpha = 0.3) +
    labs(y = "Predicted Viral Sharing", x = "Phylogenetic Similarity") +
    scale_colour_continuous_sequential(palette = AlberPalettes[1]) +
    theme(legend.position = "top") +
    geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
  
  FitList[["VirusBinary"]] %>% filter(!Phylo == last(unique(FitList[["VirusBinary"]]$Phylo))&
                                        DietSim == last(unique(FitList[["VirusBinary"]]$DietSim))) %>%
    
    ggplot(aes(Space, Fit, colour = Phylo)) + 
    geom_line(aes(group = as.factor(Phylo)), alpha = 0.3) +
    labs(y = "Predicted Viral Sharing", x = "Geographic Overlap") +
    scale_colour_continuous_sequential(palette = AlberPalettes[2]) +
    theme(legend.position = "top") +
    geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01)
  
) %>% arrange_ggplot2

dev.off()

# Figure 2.	Observed hosts have higher Predicted degree centrality #####

Errordf <- Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% group_by(Obs) %>%
  summarise(CentreDegree = mean(ScalePredDegree),
            sd = sd(ScalePredDegree),
            N = n()) %>% mutate(se = sd/(N^0.5))

Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% 
  ggplot(aes(as.factor(Subset), ScalePredDegree, colour = as.factor(Subset))) + 
  #geom_violin() +
  geom_sina(aes(alpha = as.factor(Subset))) + scale_alpha_manual(values = c(0.3,1,1,1)) +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  theme(legend.position = "none") +
  geom_point(data = Errordf, colour = "black", aes(y = CentreDegree)) + 
  geom_errorbar(data = Errordf, inherit.aes = F, 
                aes(x = as.factor(Subset), 
                    ymin = CentreDegree - se,
                    ymax = CentreDegree + se), width = 0.1) +
  scale_x_discrete(labels = c("Unobserved", "HP3", "EID", "Both")) +
  scale_colour_manual(values = c("grey", AlberColours[[2]], AlberColours[[1]], AlberColours[[3]])) +
  ggsave("Figures/Order_ObsYN_Scale_Sina.jpeg", units = "mm", height = 150, width = 150, dpi = 600)

# Figure 3.	Mammal order level centrality #####

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


# Figure 4.	Gridcell level centrality measure #####

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

ggplot(Hosts, aes(Degree, PredDegree1, colour = Sp)) + 
  geom_point(alpha = 0.5) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(x = "Observed Degree", y = "Predicted Degree", title = "With Random Effects") +
  lims(x = c(0, max(Hosts$Degree, na.rm = T)), y = c(0, max(Hosts$Degree, na.rm = T))) +
  scale_colour_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("SIFigures/Degree.PredDegree1.jpeg", units = "mm", width = 100, height = 100)

ggplot(Hosts, aes(Degree, PredDegree1b, colour = Sp)) + 
  geom_point(alpha = 0.5) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(x = "Observed Degree", y = "Predicted Degree", title = "No Random Effects") +
  lims(x = c(0, max(Hosts$Degree, na.rm = T)), y = c(0, max(Hosts$Degree, na.rm = T))) +
  scale_colour_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("SIFigures/Degree.PredDegree1b.jpeg", units = "mm", width = 100, height = 100)

jpeg("SIFigures/DegreePairs.jpeg", units = "mm", width = 200, height = 200, res = 300)

GGally::ggpairs(Hosts %>% select(contains("Degree")), 
                lower = list(continuous = wrap("smooth", colour = AlberColours[2])))

dev.off()

# Order-level Correlations among degree predictions ####

BarBarGraph(Hosts, "hOrder", "Degree", "AllPredDegree") +
  ggtitle("Order-Level Predictions Do Not Correlate") +
  ggsave("SIFigures/OrderDegreePredicts.jpeg", units = "mm", height = 150, width = 150, dpi = 300)

# Species in the observed dataset have higher predicted centrality across all mammals ####

Panth1 %>% 
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>%
  BarGraph(., "hOrder", "AllPredDegree", "Obs", text = "N", order = T) +
  labs(fill = "Observed") +
  ggtitle("Species in the HP3 dataset have higher predicted centrality across all mammals") +
  #scale_x_discrete(limits = ) +
  ggsave("SIFigures/Order_ObsYN.jpeg", units = "mm", height = 100, width = 200)

Panth1 %>% 
  filter(hOrder %in% (Panth1 %>% filter(EIDObs==1) %>% droplevels)$hOrder) %>%
  BarGraph(., "hOrder", "AllPredDegree", "EIDObs", text = "N", order = T) +
  labs(fill = "Observed") +
  ggtitle("Species in the EID dataset have higher predicted centrality across all mammals") +
  #scale_x_discrete(limits = ) +
  ggsave("SIFigures/Order_EIDObsYN.jpeg", units = "mm", height = 100, width = 200)

Panth1 %>% 
  filter(hOrder %in% (Panth1 %>% filter(EIDObs==1) %>% droplevels)$hOrder) %>%
  ggplot(aes(hOrder, AllPredDegree, colour = as.factor(EIDObs))) + 
  ggforce::geom_sina(position = position_dodge(w = 0.5))

Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% group_by(Obs) %>%
  summarise(CentreDegree = mean(ScalePredDegree),
            sd = sd(ScalePredDegree),
            N = n()) %>% mutate(se = sd/(N^0.5)) %>% ggplot(aes(as.factor(Obs), CentreDegree, fill = as.factor(Obs))) + geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = CentreDegree-se, ymax = CentreDegree+se), width = 0.1)  +
  theme(legend.position = "none") +
  labs(x = "Observed Data", y = "Within-order scaled degree (SD)") +
  ggsave("Figures/Order_ObsYN_Scale.jpeg", units = "mm", height = 100, width = 100)

# Bats specifically ####

Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder == "Chiroptera") %>% droplevels %>%
  ggplot(aes(as.factor(Subset), ScalePredDegree, colour = as.factor(Subset))) + 
  #geom_violin() +
  geom_sina() +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  ggtitle("Non-central bats are poorly sampled") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Unobserved", "HP3", "EID", "Both")) +
  scale_colour_manual(values = c("grey", AlberColours[[2]], AlberColours[[1]], AlberColours[[3]])) +
  ggsave("SIFigures/ObservedBatCentrality.jpeg", units = "mm", height = 100, width = 150)

Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder == "Chiroptera") %>% droplevels %>% 
  filter(!MSW05_Family %in% c("Craseonycteridae", "Furipteridae", "Noctilionidae", "Thyropteridae")) %>%
  ggplot(aes(as.factor(Subset), ScalePredDegree, colour = as.factor(Subset))) + 
  geom_sina() +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  theme(legend.position = "none") +
  ggtitle("Non-central families of bats are poorly sampled") +
  scale_x_discrete(labels = c("Unobserved", "HP3", "EID", "Both")) +
  scale_colour_manual(values = c("grey", AlberColours[[2]], AlberColours[[1]], AlberColours[[3]])) +
  facet_wrap(~MSW05_Family) +
  ggsave("SIFigures/ObservedBatFamilyCentrality.jpeg", units = "mm", height = 100, width = 200)

# Rodents Specifically ####

Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder == "Rodentia") %>% droplevels %>%
  ggplot(aes(as.factor(Subset), ScalePredDegree, colour = as.factor(Subset))) + 
  #geom_violin() +
  geom_sina() +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  ggtitle("Rodents are evenly sampled") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Unobserved", "HP3", "EID", "Both")) +
  scale_colour_manual(values = c("grey", AlberColours[[2]], AlberColours[[1]], AlberColours[[3]])) +
  ggsave("SIFigures/ObservedRodentCentrality.jpeg", units = "mm", height = 100, width = 150)

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

OrderPairs %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1), colour = MisspentEuth)) + 
  coord_fixed() +
  geom_point(alpha = 0.6) +
  labs(x = "log(Order 1 Hosts) + log(Order 2 Hosts)", 
       y = "log(Predicted Links + 1)",
       title = "Scaling of between-order links",
       colour = "Eutherian Orders") +
  scale_colour_manual(values = c(AlberColours[[3]], AlberColours[[2]],AlberColours[[1]])) +
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

# No other viral traits matter for prediction ####

jpeg("SIFigures/ViralTraits_Predictability.jpeg", units = "mm", width = 200, height = 200, res = 300)

VirusCovar %>% 
  lapply(function(a){
    BarGraph(ValidSummary, a, "MeanRank", text = "N") +
      theme(legend.position = "none")
  }) %>% 
  arrange_ggplot2(ncol = 3)

dev.off()



