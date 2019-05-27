
# Generating Figures ####

library(tidyverse); library(ggregplot); library(colorspace); library(cowplot)

# Figure 1: GAMM model outputs ####

load("Output Files/FitList.Rdata")
load("Output Files/BAMList.Rdata")

plot_grid( FitList[["VirusBinary"]] %>% 
             filter(!is.na(SpaceQuantile)) %>%
             ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + 
             geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
             geom_line(aes(group = as.factor(Space))) +
             labs(y = "Viral Sharing Probability", x = "Phylogenetic Similarity", 
                  colour = "Overlap", fill = "Overlap") +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
             scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
             theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
             geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
           
           FitList[["VirusBinary"]] %>% 
             filter(!is.na(PhyloQuantile)) %>%
             ggplot(aes(Space, Fit, colour = PhyloQuantile)) + 
             geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
             geom_line(aes(group = as.factor(Phylo))) +
             labs(y = "Viral Sharing Probability", x = "Geographic Overlap", 
                  colour = "Relatedness", fill = "Relatedness") +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
             scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
             theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
             geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01),
           
           FitList[["VirusBinary"]] %>% 
             filter(!Phylo == last(unique(Phylo)),
                    !Space == last(unique(Space))) %>%
             ggplot(aes(Space, Phylo)) + 
             geom_tile(aes(fill = Fit)) + 
             labs(x = "Geographic Overlap", 
                  y = "Phylogenetic Similarity",
                  fill = "Estimate") +
             #ggtitle("Tensor Field") +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             theme(legend.position = "bottom") +
             scale_fill_continuous_sequential(palette = "Greens 2", cmax = 20, end = 1,
                                              limits = c(0,1),
                                              breaks = c(0,0.5,1)),
           
           DataList$VirusBinary %>%
             ggplot(aes(Space, Phylo)) + 
             labs(x = "Geographic Overlap", 
                  y = "Phylogenetic Similarity") +
             #ggtitle("Data Distribution") +
             scale_fill_continuous_sequential(palette = "purp", begin = 0.2) +
             lims(x = c(0,1), y = c(0,1)) +
             coord_fixed() +
             theme(legend.position = "bottom") +
             geom_hex(aes(fill = stat(log(count)))),
           
           nrow = 2, 
           rel_heights = c(1,1.23), 
           labels = "AUTO") %>% 
  save_plot(filename = "Figures/Figure1.jpeg", 
            #units = "mm", width = 200, height = 200,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1)

# Figure 2.	Observed hosts have higher predicted degree centrality #####

Errordf <- Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% 
  group_by(Subset) %>%
  summarise(CentreDegree = mean(ScalePredDegree, na.rm = T),
            sd = sd(ScalePredDegree, na.rm = T),
            N = n()) %>% mutate(se = sd/(N^0.5))

Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% 
  ggplot(aes(Subset, ScalePredDegree, colour = as.factor(Subset))) + 
  # geom_violin(colour = "grey") +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_sina(aes(alpha = Subset)) + scale_alpha_manual(values = c(0.3,1,1,1)) +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  theme(legend.position = "none") +
  geom_point(data = Errordf, colour = "black", aes(y = CentreDegree)) + 
  geom_errorbar(data = Errordf, inherit.aes = F, 
                aes(x = as.factor(Subset), 
                    ymin = CentreDegree - se,
                    ymax = CentreDegree + se), width = 0.1) +
  scale_x_discrete(labels = c("Unobserved", "EID Only","HP3 Only",  "Both")) +
  scale_colour_manual(values = c("light grey", AlberColours[[1]], AlberColours[[2]], AlberColours[[3]])) +
  ggsave("Figures/Order_ObsYN_Scale_Sina.jpeg", units = "mm", height = 150, width = 150, dpi = 600)

plot1 <- Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% 
  ggplot(aes(Subset, ScalePredDegree, colour = as.factor(Subset))) + 
  # geom_violin(colour = "grey") +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_sina(aes(alpha = Subset)) + scale_alpha_manual(values = c(0.3,1,1,1)) +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  theme(legend.position = "none") +
  geom_point(data = Errordf, colour = "black", aes(y = CentreDegree)) + 
  geom_errorbar(data = Errordf, inherit.aes = F, 
                aes(x = as.factor(Subset), 
                    ymin = CentreDegree - se,
                    ymax = CentreDegree + se), width = 0.1) +
  scale_x_discrete(labels = c("Unobserved", "EID Only","HP3 Only",  "Both")) +
  scale_colour_manual(values = c("light grey", AlberColours[[1]], AlberColours[[2]], AlberColours[[3]]))


Hosts$AnyZoo <- as.factor(as.numeric(Hosts$hZoonosisCount>0))

Hosts <- Hosts %>% left_join(AllDegrees)

Hosts %>% # filter(AllPredDegree<250) %>% 
  BarGraph("AnyZoo", "AllPredDegree", Text = "N") + 
  labs(x = "Zoonotic Host", y = "Predicted links") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c(AlberColours[[1]],AlberColours[[2]])) +
  ggsave("Figures/Zoonosis_PredictedLinks.jpeg", units = "mm", width = 100, height = 100, dpi = 600)

plot2 <- Hosts %>% 
  BarGraph("AnyZoo", "AllPredDegree", Text = "N") + 
  labs(x = "Zoonotic Host", y = "Predicted links") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c(AlberColours[[1]],AlberColours[[2]])) 

# Species that share in EID have higher sharing probability in our networks ####

SinaGraph(EIDCordf, "EIDConnected", "PredNetwork", Scale = "width", Alpha = 0.2) +
  #scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[3]], nmax = 8, order = c(4,7)) +
  scale_alpha_manual(values = c(0.2,0.2)) +
  theme(legend.position = "none") +
  labs(x = "EID2 Sharing", y = "Predicted Sharing Probability") +
  ggsave("Figures/EIDConnections_Predictions.jpeg", units = "mm", height = 100, width = 150, dpi = 300)

table(EIDCordf$EIDConnected)

plot3 <-
  SinaGraph(EIDCordf, "EIDConnected", "PredNetwork", Scale = "width", Alpha = 0.2) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  # scale_colour_discrete_sequential(palette = AlberPalettes[[3]], nmax = 8, order = c(4,7)) +
  scale_alpha_manual(values = c(0.2,0.2)) +
  theme(legend.position = "none") +
  labs(x = "EID2 Sharing", y = "Predicted Sharing Probability") +
  lims(y = c(0,1)) + geom_text(data = data.frame(),
                               inherit.aes = F, 
                               aes(label = c("N=46887", "N=7069"), 
                                   x = as.factor(c(0,1)), y = c(1,1)))

# Puttin these verification figures together in one ####

bottom_row <- plot_grid(plot2, plot3, labels = c("B","C"))

plot_grid(plot1, bottom_row, nrow = 2, 
          labels = c("A",NULL), 
          rel_heights = c(1.5,1)) %>%
  save_plot(filename = "Figures/Figure2.jpeg", 
            #units = "mm", width = 200, height = 200,
            # ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 2)

# Figure 3.	Mammal order level centrality #####

load("Output Files/Panth1.Rdata")

BarGraph(Panth1, "hOrder", "AllPredDegree", Just = T, Order = T, Text = "N") +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
  labs(x = "Host Order", y = "Predicted links", title = "All Links") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("Figures/AllPredDegree.jpeg", units = "mm", height = 120, width = 150)

BarGraph(Panth1, "hOrder", "InDegree", Just = T, Order = T, Text = "N") +
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]]) +
  labs(x = "Host Order", y = "Predicted links", title = "Within-order links") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("Figures/InDegree.jpeg", units = "mm", height = 120, width = 150)

BarGraph(Panth1, "hOrder", "OutDegree", Just = T, Order = T, Text = "N") +
  scale_fill_discrete_sequential(palette = AlberPalettes[[3]]) +
  labs(x = "Host Order", y = "Predicted links", title = "Out-of-order links") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("Figures/OutDegree.jpeg", units = "mm", height = 120, width = 150)

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
  coord_fixed() + 
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  theme_void() +
  ggsave("Figures/All Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

PlotGrids %>% filter(Metric == "InDegree")  %>% mutate(Degree = ifelse(Degree>200, 200, ifelse(Degree<30,40,Degree))) %>%
  ggplot(aes(x, y, fill = Degree)) + # , colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  coord_fixed() +  
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  theme_void() +
  ggsave("Figures/In Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

PlotGrids %>% filter(Metric == "OutDegree")  %>% mutate(Degree = ifelse(Degree>150, 150, ifelse(Degree<110,110,Degree))) %>%
  ggplot(aes(x, y, fill = Degree)) + #, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  coord_fixed() +  
  scale_colour_continuous_sequential(palette = AlberPalettes[3]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[3]) +
  theme_void() +
  ggsave("Figures/Out Link Map.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

# Putting together geog and taxon predictions ####

TaxonPlot <- BarGraph(Panth1, "hOrder", "AllPredDegree", Just = T, Order = T, Text = "N") +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
  labs(x = "Host Order", y = "Predicted links", title = "All Links") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

GeogPlot <- PlotGrids %>% filter(Metric == "AllPredDegree") %>% mutate(Degree = ifelse(Degree>320, 320, Degree)) %>%
  ggplot(aes(x, y, fill = Degree)) + #, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  coord_fixed() + 
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  theme_void()

plot_grid(TaxonPlot + labs(x = NULL), 
          GeogPlot, 
          nrow = 2, 
          rel_heights = c(1, 1), 
          labels = "AUTO") %>%
  save_plot(filename = "Figures/Figure2.jpeg", 
            #units = "mm", width = 200, height = 200,
            # ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 2)

# Trying just one ####

Map_All <- PlotGrids %>% filter(Metric == "AllPredDegree") %>% mutate(Degree = ifelse(Degree>320, 320, Degree)) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) +
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  coord_fixed() + lims(x = c(225, 1350)) + 
  guides(alpha = "none") + 
  labs(fill = "All links", colour = "All links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  theme_void() + 
  theme(legend.position = "bottom")

Map_In <- PlotGrids %>% filter(Metric == "InDegree")  %>% mutate(Degree = ifelse(Degree>200, 200, ifelse(Degree<30,40,Degree))) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  coord_fixed() + lims(x = c(225, 1350)) +  
  guides(alpha = "none")+ 
  labs(fill = "Within-order links", colour = "Within-order links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  theme_void() + 
  theme(legend.position = "bottom") 

Map_Out <- PlotGrids %>% filter(Metric == "OutDegree")  %>% 
  mutate(Degree = ifelse(Degree>150, 150, ifelse(Degree<110,110,Degree))) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log(Density))) +
  coord_fixed() + lims(x = c(225, 1350)) +  
  guides(alpha = "none")+ 
  labs(fill = "Out-of-order links", colour = "Out-of-order links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[3]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[3]) +
  theme_void() + 
  theme(legend.position = "bottom") 

TextSize = 3
AxisTextX = 8
AxisTextY = 10

Taxon_All <- BarGraph(Panth1, "hOrder", "AllPredDegree", Just = T, Order = T, Text = "N", TextSize = TextSize) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
  labs(x = NULL, y = "All links") +
  theme(axis.text.x = element_text(size = AxisTextX),
        axis.text.y = element_text(size = AxisTextY)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

Taxon_In <- BarGraph(Panth1, "hOrder", "InDegree", Just = T, Order = T, Text = "N", TextSize = TextSize) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]]) +
  labs(x = NULL, y = "Within-order links") +
  theme(axis.text.x = element_text(size = AxisTextX),
        axis.text.y = element_text(size = AxisTextY)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

Taxon_Out <- BarGraph(Panth1, "hOrder", "OutDegree", Just = T, Order = T, Text = "N", TextSize = TextSize) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[3]]) +
  labs(x = NULL, y = "Out-of-order links") +
  theme(axis.text.x = element_text(size = AxisTextX),
        axis.text.y = element_text(size = AxisTextY)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

Row_All <- plot_grid(Taxon_All, Map_All, 
                     nrow = 1, rel_widths = c(1.5,1), 
                     labels = LETTERS[1:2], axis = "t")

Row_In <- plot_grid(Taxon_In, Map_In, nrow = 1, rel_widths = c(1.5,1), 
                    labels = LETTERS[1:2+2], axis = "t")

Row_Out<- plot_grid(Taxon_Out, Map_Out, nrow = 1, rel_widths = c(1.5,1), 
                    labels = LETTERS[1:2+4], axis = "t")

WholePlot <- plot_grid(Row_All, Row_In, Row_Out, nrow = 3)

WholePlot %>%  save_plot(filename = "Figures/Figure3.jpeg", 
                         #units = "mm", width = 200, height = 200,
                         # ncol = 2, # we're saving a grid plot of 2 columns
                         # nrow = 2, # and 2 rows
                         # each individual subplot should have an aspect ratio of 1.3
                         base_aspect_ratio = 1,
                         base_height = 9)

# 6.	Panel of maps: Illustrative maps of host ranges for mammal species from host predictions (currently unidentified) for specific viruses, e.g. CCHFV; Nipah virus; Ebola; etc. (supplement will list out probability of being a host for all viruses)

PredPlot(HostList = VirusAssocs[["Crimean-Congo_hemorrhagic_fever_virus"]], 
         Focal = c("Observed"), 
         Facet = F,
         Validate = F)$MapPlot + 
  theme_void() + theme(plot.title = element_text(hjust=0.5)) +
  ggtitle("CCHF Known Hosts") +
  ggsave("Figures/CCHF_Known.jpeg", units = "mm", width = 250, height = 150)

PredPlot(HostList = VirusAssocs[["Crimean-Congo_hemorrhagic_fever_virus"]], 
         Focal = c("Predicted"), 
         Facet = F,
         Validate = F, Threshold = 15)$MapPlot + 
  theme_void() + theme(plot.title = element_text(hjust=0.5)) +
  #ggtitle("CCHF Model Success") +
  ggsave("Figures/CCHF_Predicted.jpeg", units = "mm", width = 250, height = 150)

PredPlot(HostList = VirusAssocs[["Crimean-Congo_hemorrhagic_fever_virus"]], 
         Focal = c("Predicted","Observed"), 
         Facet = T,
         Validate = F, Summarise = F)$MapPlot + 
  ggtitle("CCHF Model Success") +
  ggsave("Figures/CCHF_Success.jpeg", units = "mm", width = 100, height = 100)

PredPlot(Virus = "Crimean-Congo_hemorrhagic_fever_virus", 
         Focal = c(0), 
         Validate = F,
         Facet = F)

# Supplementary figures #####

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

# Subnetwork model outputs ####

# Phylo ####

RespLabels <- c("Viruses", "RNA", "DNA", "Vector-borne", "Non-vector")

lapply(2:5, function(a){
  
  FitList[[a]] %>% 
    filter(!is.na(SpaceQuantile)) %>%
    ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
    geom_line(aes(group = as.factor(Space))) +
    labs(y = "Predicted Viral Sharing", x = "Phylogenetic Similarity", 
         colour = "Overlap", fill = "Overlap",
         title = RespLabels[a]) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
    theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
    geom_rug(data = DataList[[a]], inherit.aes = F, aes(x = Phylo), alpha = 0.01)
  
  
}) %>% plot_grid(plotlist = .) %>%   
  save_plot(filename = "SIFigures/SubModels_Phylo.jpeg", 
            #units = "mm", width = 200, height = 200,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1)


# Space ####

lapply(2:5, function(a){
  
  FitList[[a]] %>% 
    filter(!is.na(PhyloQuantile)) %>%
    ggplot(aes(Space, Fit, colour = PhyloQuantile)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
    geom_line(aes(group = as.factor(Phylo))) +
    labs(y = "Predicted Viral Sharing", x = "Geographic Overlap", 
         colour = "Relatedness", fill = "Relatedness",
         title = RespLabels[a]) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
    scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
    theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
    geom_rug(data = DataList[[a]], inherit.aes = F, aes(x = Space), alpha = 0.01)
  
}) %>% plot_grid(plotlist = .) %>%   
  save_plot(filename = "SIFigures/SubModels_Space.jpeg", 
            #units = "mm", width = 200, height = 200,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1)



# Tensor 

lapply(2:5, function(a){
  
  FitList[[a]] %>% 
    filter(!Phylo == last(unique(Phylo)),
           !Space == last(unique(Space))) %>%
    ggplot(aes(Space, Phylo)) + 
    geom_tile(aes(fill = Fit)) + 
    labs(x = "Geographic Overlap", 
         y = "Phylogenetic Similarity",
         fill = "Estimate",
         title = Resps[a]) +
    #ggtitle("Tensor Field") +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    theme(legend.position = "bottom") +
    scale_fill_continuous_sequential(palette = "Greens 2", cmax = 20, end = 1)
  
}) %>% plot_grid(plotlist = .) %>%   
  save_plot(filename = "SIFigures/SubModels_Tensor.jpeg", 
            #units = "mm", width = 200, height = 200,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1)



lapply(2:5, function(a){
  
  DataList[[a]] %>% 
    ggplot(aes(Space, Phylo)) + 
    labs(x = "Geographic Overlap", 
         y = "Phylogenetic Similarity",
         title = Resps[a]) +
    #ggtitle("Data Distribution") +
    scale_fill_continuous_sequential(palette = "purp", begin = 0.2) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = stat(log(count))))
  
}) %>% plot_grid(plotlist = .) %>%   
  save_plot(filename = "SIFigures/SubModels_Distributions.jpeg", 
            #units = "mm", width = 200, height = 200,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1)


# Data Description: observed tiles of space:phylo with virusbinary fill #####

FinalHostMatrix %>%
  mutate(SpaceQuantile = cut(Space, breaks = c(0:10)/10, include.lowest = T),
         PhyloQuantile = cut(Phylo, breaks = c(0:10)/10, include.lowest = T)) %>% 
  group_by(SpaceQuantile, PhyloQuantile) %>% summarise(VirusBinary = mean(VirusBinary), N = n()) %>% 
  ggplot(aes(SpaceQuantile, PhyloQuantile, fill = VirusBinary, alpha = log10(N))) + 
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]]) +
  coord_fixed() +
  geom_tile() + geom_text(aes(label = N), size = 2.5)

# Data Description: observed versus all-mammal data distributions #####

plot_grid(
  
  DataList$VirusBinary %>%
    ggplot(aes(Space, Phylo)) + 
    labs(x = "Geographic Overlap", 
         y = "Phylogenetic Similarity") +
    #ggtitle("Data Distribution") +
    scale_fill_continuous_sequential(palette = "purp", begin = 0.2) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = stat(log(count)))),
  
  AllMammaldf %>%
    ggplot(aes(Space, Phylo)) + 
    labs(x = "Geographic Overlap", 
         y = "Phylogenetic Similarity") +
    #ggtitle("Data Distribution") +
    scale_fill_continuous_sequential(palette = "purp", begin = 0.2) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = stat(log(count)))),
  
  ncol = 2, 
  rel_heights = c(1,1), 
  labels = "AUTO") %>% 
  save_plot(filename = "SIFigures/Data Distributions.jpeg", 
            ncol = 2, # we're saving a grid plot of 2 columns
            base_aspect_ratio = 1)

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
  mutate(Obs = as.factor(Obs)) %>%
  # filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>%
  BarGraph(., "hOrder", "AllPredDegree", "Obs", Text = "N", Order = T, Just = T) +
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
            N = n()) %>% mutate(se = sd/(N^0.5)) %>% 
  ggplot(aes(as.factor(Obs), CentreDegree, fill = as.factor(Obs))) + 
  geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = CentreDegree-se, ymax = CentreDegree+se), width = 0.1)  +
  theme(legend.position = "none") +
  labs(x = "Observed Data", y = "Within-order scaled links") +
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

# Scaling of within- and between-order networks ####

Panth1 %>% group_by(hOrder) %>%
  summarise(Number = n(),
            AllPredDegree = mean(AllPredDegree),
            InDegree = mean(InDegree),
            OutDegree = mean(OutDegree)) %>% lm(log(InDegree+1) ~ log(Number+1), data = .)

OrderLevelLinks %>% ggplot(aes(log(HostNumber), log(Degree+1))) + 
  geom_point(colour = AlberColours[[3]]) + 
  coord_fixed() + 
  geom_smooth(colour = "black", fill = NA, method = lm) + 
  stat_smooth(fill = NA, geom = "ribbon", lty = 2, colour = "black", method = lm) +
  facet_wrap(~Metric, labeller = labeller(Metric = c("AllPredDegree" = "All links",
                                                     "OutDegree" = "Out-of-order links",
                                                     "InDegree" = "Within-order links")))  +
  theme(strip.background = element_rect(fill = "white", colour = "grey")) +
  labs(x = "log(Host Number)") + 
  ggsave("SIFigures/log_NHosts_Degree_Order.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

OrderPairs %>%
  ggplot(aes(log(`1`) + log(`2`), log(Degree+1))) + 
  coord_fixed() +
  geom_point(alpha = 0.6) +
  labs(x = "log(order 1 hosts*order 2 hosts)", 
       y = "log(Predicted links + 1)",
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

# Correlations among host range and predictability ####

ValidSummary %>% 
  ggplot(aes(HostRangeMean, log10(MeanRank))) +
  geom_point(aes(colour = vFamily)) + 
  geom_smooth(colour = "black", fill = NA, method = lm) +
  stat_smooth(fill = NA, geom = "ribbon", lty = 2, colour = "black", method = lm) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[2]]) + 
  labs(x = "Average host phylogenetic similarity",
       y = "log10(mean rank)") +
  ggsave("SIFigures/HostRange_Predictability.jpeg", units = "mm", 
         width = 200, height = 100, dpi = 300)

# RNA Viruses are harder to predict ####

BarGraph(ValidSummary, "vDNAoRNA", "MeanRank", Text = "N") +
  scale_fill_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  theme(legend.position = "none") +
  ggsave("SIFigures/RNA_DNA_Prediction.jpeg", units = "mm", height = 100, width = 100)

# No other viral traits matter for prediction ####

jpeg("SIFigures/ViralTraits_Predictability.jpeg", units = "mm", width = 200, height = 200, res = 300)

VirusCovar %>% 
  lapply(function(a){
    BarGraph(ValidSummary, a, "MeanRank", Text = "N") +
      theme(legend.position = "none")
  }) %>% 
  arrange_ggplot2(ncol = 3)

dev.off()

# Taxonomic patterns of predictability ####

Errordf <- ValidSummary %>% group_by(vFamily) %>% 
  summarise(MedianRank = median(MeanRank),
            sd = sd(MeanRank),
            N = n()) %>% mutate(se = sd/(N^0.5)) %>%
  slice(order(MedianRank))

vFamilyOrder <- Errordf$vFamily

ggplot(ValidSummary, aes(vFamily, log10(MeanRank))) + geom_sina(colour = AlberColours[[2]]) + 
  scale_x_discrete(limits = vFamilyOrder) +
  labs(x = "Viral Family", y = "log10(Focal host rank)", colour = "Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("SIFigures/VirusTaxonomy_PredictionSuccess.jpeg", units = "mm", width = 150, height = 100, dpi = 300)

# Simple Presentation stuff ####

ggtree(chiroptera, aes(color = group, alpha = group)) +
  scale_colour_manual(values = c("black", "red", "blue")) +
  scale_alpha_manual(values = c(0.01,0.5,1)) +
  ggsave("Tree.jpeg", units = "mm", width = 100, height = 200, dpi = 300)


# Correlation between mean rank of focal host predictions and the proportion of links they're present for

ggplot(ValidSummary, aes(log10(MeanRank), MeanCount1)) + geom_point() + 
  geom_smooth(method = lm) +
  ggsave("SIFigures/Rank_Count.jpeg", units = "mm", height = 100, width = 100)

BarGraph(ValidSummary, "vFamily", "MeanRank", Text = "N", Order = T, Just = T) +
  labs(x = "Viral Family")

BarGraph(ValidSummary, "vFamily", "LogRank", Text = "N", Order = T, Just = T) +
  labs(x = "Viral Family")





