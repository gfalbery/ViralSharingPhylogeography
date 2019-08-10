
# Generating Figures ####

library(tidyverse); library(ggregplot); library(colorspace); library(cowplot)

theme_set(theme_cowplot())

# Figure 1: GAMM model outputs ####

load("Output Files/FitList.Rdata")
load("Output Files/BAMList.Rdata")

plot_grid(FitList[["VirusBinary"]] %>% 
            filter(!is.na(SpaceQuantile)) %>%
            ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
            geom_line(aes(group = as.factor(Space))) +
            labs(y = "Viral sharing probability", x = "Phylogenetic similarity", 
                 colour = "Geographic overlap", fill = "Geographic overlap") +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
            scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
            theme(legend.position = c(0.1, 0.8), 
                  legend.title = element_text(size = 10),
                  legend.background = element_rect(colour = "dark grey")) +
            geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
          
          FitList[["VirusBinary"]] %>% 
            filter(!is.na(PhyloQuantile)) %>%
            ggplot(aes(Space, Fit, colour = PhyloQuantile)) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
            geom_line(aes(group = as.factor(Phylo))) +
            labs(y = "Viral sharing probability", x = "Geographic overlap", 
                 colour = "Phylogenetic similarity", fill = "Phylogenetic similarity") +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
            scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
            theme(legend.position = c(0.1, 0.8), 
                  legend.title = element_text(size = 10),
                  legend.background = element_rect(colour = "dark grey")) +
            geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01),
          
          FitList[["VirusBinary"]] %>% 
            filter(!Phylo == last(unique(Phylo)),
                   !Space == last(unique(Space))) %>%
            ggplot(aes(Space, Phylo)) + 
            geom_tile(aes(fill = Fit)) + 
            labs(x = "Geographic overlap", 
                 y = "Phylogenetic similarity",
                 fill = "Viral sharing\nprobability") +
            #ggtitle("Tensor Field") +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            theme(legend.position = "bottom",
                  legend.title = element_text(size = 10)) +
            geom_contour(aes(z = Fit), colour = "white", alpha = 0.8) + 
            metR::geom_text_contour(aes(z = Fit), colour = "white", size = 2.5, hjust = 0.5, vjust = 1.1, check_overlap = T) +
            scale_fill_continuous_sequential(palette = "ag_GrnYl",
                                             limits = c(0,1),
                                             breaks = c(0,0.5,1)),
          
          DataList$VirusBinary %>%
            ggplot(aes(Space, Phylo)) + 
            labs(x = "Geographic overlap", 
                 y = "Phylogenetic similarity") +
            #ggtitle("Data Distribution") +
            scale_fill_continuous_sequential(palette = "Heat 2") +
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

load("Output Files/Panth1.Rdata")

Errordf <- Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% 
  group_by(Subset) %>%
  summarise(CentreDegree = mean(ScalePredDegree, na.rm = T),
            sd = sd(ScalePredDegree, na.rm = T),
            N = n()) %>% mutate(se = sd/(N^0.5))

plot1 <- Panth1 %>% group_by(hOrder) %>%
  mutate(ScalePredDegree = scale_this(AllPredDegree)) %>%   
  filter(hOrder %in% (Panth1 %>% filter(Obs==1) %>% droplevels)$hOrder) %>% 
  ggplot(aes(Subset, ScalePredDegree, colour = as.factor(Subset))) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_sina(aes(alpha = Subset)) + scale_alpha_manual(values = c(0.3,0.7,0.7,0.7)) +
  labs(x = "Host dataset", y = "Within-order scaled degree (SD)") +
  theme(legend.position = "none") +
  geom_point(data = Errordf, colour = "black", aes(y = CentreDegree)) + 
  geom_errorbar(data = Errordf, inherit.aes = F, 
                aes(x = as.factor(Subset), 
                    ymin = CentreDegree - se,
                    ymax = CentreDegree + se), 
                width = 0.1) +
  scale_x_discrete(labels = c("Unobserved", "EID2 only","Training data only",  "Both")) +
  scale_colour_manual(values = c("light grey", AlberColours[[1]], AlberColours[[2]], AlberColours[[3]]))

Hosts$AnyZoo <- as.factor(as.numeric(Hosts$hZoonosisCount>0))

Ns1 = Hosts %>% filter(Sp%in%FHN) %>% pull(AnyZoo) %>% table()

plot2 <- Hosts %>% 
  SinaGraph("AnyZoo", 
            "AllPredDegree.x",
            Alpha = 0.6) + 
  labs(x = "Zoonotic host", 
       y = "Predicted links") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  scale_x_discrete(labels = c("No","Yes")) +
  geom_text(data = data.frame(),
            inherit.aes = F, 
            aes(label = paste0("N=",c(Ns1[1],Ns1[2])),
                x = as.factor(c(0,1)), y = rep(575,2)))

EIDCordf %>% 
  filter(!(Sp%in%FHN&Sp2%in%FHN)) -> 
  SubEIDCordf

Ns = table(SubEIDCordf$EIDConnected)

plot3 <-
  SinaGraph(SubEIDCordf, 
            "EIDConnected", 
            "PredNetwork", 
            Scale = "width", Alpha = 0.2) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  scale_alpha_manual(values = c(0.2,0.2)) +
  scale_x_discrete(labels = c("No", "Yes")) +
  theme(legend.position = "none") +
  labs(x = "Shares in EID2", 
       y = "Sharing probability") +
  lims(y = c(0,1)) + 
  geom_text(data = data.frame(),
            inherit.aes = F, 
            aes(label = paste0("N=",c(Ns[1],Ns[2])), 
                x = as.factor(c(0,1)), y = c(1,1)))

bottom_row <- plot_grid(plot2, plot3, labels = c("B","C"))

plot_grid(plot1, bottom_row, nrow = 2, 
          labels = c("A",NULL), 
          rel_heights = c(1.5,1)) %>%
  save_plot(filename = "Figures/Figure2a.jpeg", 
            nrow = 2, # and 2 rows
            base_aspect_ratio = 2)

bottom_row <- plot_grid(plot3, plot2, labels = c("A","B"))

plot_grid(bottom_row, plot1, nrow = 2, 
          labels = c(" ","C"), 
          rel_heights = c(1,1.5)) %>%
  save_plot(filename = "Figures/Figure2b.jpeg", 
            nrow = 2, # and 2 rows
            base_aspect_ratio = 2)

# Figure 3.	taxonomic and geographic prediction patterns #####

load("Output Files/Panth1.Rdata")
load("Output Files/GridDegree.Rdata")

Map_All <- GridDegree %>% filter(Metric == "AllDegree") %>% 
  mutate(Degree = ifelse(Degree>300, 300, Degree)) %>%
  mutate(RichCut = as.factor(ifelse(Richness>2,1,0))) %>%
  ggplot(aes(X, Y, fill = Degree, colour = Degree)) +
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log10(Richness+1))) +
  coord_fixed() + 
  guides(alpha = "none") +
  labs(fill = "All links", colour = "All links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  theme_void() + 
  theme(legend.position = "bottom")

Map_In <- GridDegree %>% filter(Metric == "InDegree")  %>% 
  mutate(Degree = ifelse(Degree>220, 220, Degree)) %>%
  mutate(RichCut = as.factor(ifelse(Richness>2,1,0))) %>%
  ggplot(aes(X, Y, fill = Degree, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log10(Richness+1))) +
  coord_fixed() + 
  guides(alpha = "none") +
  labs(fill = "Within-order links", colour = "Within-order links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  theme_void() + 
  theme(legend.position = "bottom") 

Map_Out <- GridDegree %>% filter(Metric == "OutDegree")  %>% 
  mutate(RichCut = as.factor(ifelse(Richness>2,1,0))) %>%
  mutate(Degree = ifelse(Degree>150, 150, ifelse(Degree<0, 0, Degree))) %>%
  ggplot(aes(X, Y, fill = Degree, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log10(Richness+1))) +
  coord_fixed() + 
  guides(alpha = "none") +
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
                         base_aspect_ratio = 1,
                         base_height = 9)

# Supplementary figures #####

# Species-level degree prediction correlations ####

HostsLong <- 
  Hosts %>% select(-AllPredDegree) %>% 
  gather(key = "Key", value = "Value", contains("PredDegree")) 

MaxLim <- max(HostsLong$Value, na.rm = T)

HostsLong %>%
  ggplot(aes(Degree, Value, colour = Sp)) + 
  facet_wrap(~Key, 
             labeller = labeller(Key = c("PredDegree1" = "All effects",
                                         "PredDegree1b" = "Fixed effects",
                                         "PredDegree1c" = "Random effects"))) + 
  geom_abline(lty = 2, alpha = 0.3) +
  geom_point(alpha = 0.5) +
  coord_fixed() +
  theme(legend.position = "none", strip.background = element_rect(fill = "white")) +
  labs(x = "Observed degree", y = "Predicted degree") +
  lims(x = c(0, MaxLim), y = c(0, MaxLim)) +
  scale_colour_discrete_sequential(palette = AlberPalettes[2]) +
  ggsave("SIFigures/Figure SI1.jpeg", units = "mm", width = 200, height = 100)


# Subnetwork model outputs ####

RespLabels <- c("Viruses", "RNA", "DNA", "Vector-borne", "Non-vector")

lapply(2:5, function(a){
  
  FitList[[a]] %>% 
    filter(!is.na(SpaceQuantile)) %>%
    ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
    geom_line(aes(group = as.factor(Space))) +
    labs(y = "Viral sharing probability", x = "Phylogenetic similarity", 
         colour = "Geographic overlap", fill = "Geographic overlap",
         title = RespLabels[a]) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
    theme(legend.position = c(0.1, 0.8), 
          legend.title = element_text(size = 10),
          legend.background = element_rect(colour = "dark grey")) +
    geom_rug(data = DataList[[a]], inherit.aes = F, aes(x = Phylo), alpha = 0.01)
  
  
}) %>% plot_grid(plotlist = .) %>%   
  save_plot(filename = "SIFigures/Figure SI2.jpeg", 
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            base_aspect_ratio = 1)

lapply(2:5, function(a){
  
  FitList[[a]] %>% 
    filter(!is.na(PhyloQuantile)) %>%
    ggplot(aes(Space, Fit, colour = PhyloQuantile)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
    geom_line(aes(group = as.factor(Phylo))) +
    labs(y = "Viral sharing probability", x = "Geographic overlap", 
         colour = "Phylogenetic similarity", fill = "Phylogenetic similarity",
         title = RespLabels[a]) +
    lims(x = c(0,1), y = c(0,1)) +
    coord_fixed() +
    scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
    scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
    theme(legend.position = c(0.1, 0.8), 
          legend.title = element_text(size = 10),
          legend.background = element_rect(colour = "dark grey")) +
    geom_rug(data = DataList[[a]], inherit.aes = F, aes(x = Space), alpha = 0.01)    
  
}) %>% plot_grid(plotlist = .) %>%   
  save_plot(filename = "SIFigures/Figure SI3.jpeg", 
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1)

# Species in the observed dataset have higher predicted centrality across all mammals ####

Panth1 %>% 
  filter(hOrder %in% (Panth1 %>% filter(EIDObs==1) %>% droplevels)$hOrder) %>%
  mutate(Obs = as.factor(Obs)) %>%
  BarGraph(., "hOrder", "AllPredDegree", "Obs", Text = "N", Order = T, Just = T, TextSize = 2.5) +
  labs(fill = "Observed", x  ="Host order", y = "Predicted link number") +
  scale_fill_manual(values = c(AlberColours[[1]],AlberColours[[2]]), labels = c("N", "Y")) +
  ggsave("SIFigures/Figure SI4.jpeg", units = "mm", height = 100, width = 200)

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
  labs(x = "Log(host number)", y = "Log(degree + 1)") + 
  ggsave("SIFigures/Figure SI5.jpeg", units = "mm", height = 100, width = 200, dpi = 300)

hComboList %>% group_by(Combo) %>% 
  dplyr::summarise(HostNumbers = log((c(table(Order)[1]*table(Order)[2]))),
                   Degree = log(sum(Degree)+1)) %>%
  ggplot(aes(HostNumbers, Degree)) + 
  coord_fixed() +
  geom_point(alpha = 0.6, colour = AlberColours[[3]]) +
  labs(x = "Log(order 1 hosts*order 2 hosts)", 
       y = "Log(degree + 1)") +
  ggsave("SIFigures/Figure SI6.jpeg", units = "mm", height = 100, width = 100, dpi = 300)

# Taxonomic patterns of predictability ####

ValidDF %>% group_by(vFamily) %>% 
  summarise(MedianRank = median(MeanRank),
            sd = sd(MeanRank),
            N = n()) %>% mutate(se = sd/(N^0.5)) %>%
  slice(order(MedianRank)) %>% pull(vFamily) -> vFamilyOrder

ggplot(ValidDF, aes(vFamily, log10(MeanRank))) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_boxplot() + 
  geom_point(colour = AlberColours[[2]]) + 
  scale_x_discrete(limits = vFamilyOrder) +
  labs(x = "Viral family", y = "Log10(focal host rank)", colour = "Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) + 
  scale_y_reverse() +
  ggsave("SIFigures/Figure SI7.jpeg", units = "mm", width = 150, height = 100, dpi = 300)

# Predictability is determined by host range ####

LineDF <- expand.grid(
  LogHosts = mean(ValidDF$LogHosts),
  HostRangeMean = seq(0,100)/100,
  vVectorYNna = c("Y","N"),
  vFamily = "Flaviviridae"
)

LineDF2 <- predict(LM1, 
                      newdata = LineDF,
                      se.fit = TRUE)

ValidDF %>% 
  ggplot(aes(HostRangeMean, log10(MeanRank))) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_point(aes(colour = vVectorYNna), alpha = 0.6) + 
  geom_smooth(fill = NA, method = lm, aes(colour = vVectorYNna)) +
  stat_smooth(fill = NA, aes(colour = vVectorYNna), 
              geom = "ribbon", 
              lty = 2, #colour = "black", 
              method = lm) +
  labs(x = "Average host phylogenetic similarity", y = "Log10(mean rank)", colour = "Vector-borne")  + scale_y_reverse() +
  scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[1]])) +
  ggsave("SIFigures/Figure SI8.jpeg", units = "mm", width = 150, height = 100, dpi = 300)

# Predictability is determined by host number ####

ValidDF %>% 
  ggplot(aes(LogHosts, LogRank)) +
  geom_point(alpha = 0.6, colour = AlberColours[[1]]) + 
  geom_smooth(fill = NA, method = lm, colour = "black") +
  stat_smooth(fill = NA, 
              geom = "ribbon", 
              lty = 2, colour = "black", 
              method = lm) +
  labs(x = "Log10(host number)",
       y = "Log10(mean rank)")  + #scale_y_reverse() +
  ggsave("SIFigures/HostNumber_Predictability.jpeg", units = "mm", 
         width = 150, height = 100, dpi = 300)

# No other viral traits matter for prediction ####

VirusCovar %>% 
  lapply(function(a){
    BarGraph(ValidDF, a, "MeanRank", Text = "N") +
      theme(legend.position = "none")
  }) %>% 
  arrange_ggplot2(ncol = 3)

# Summed sharing (rather than mean) ####

load("Output Files/Panth1.Rdata")
load("Output Files/GridDegree.Rdata")

Sum_Map_All <- GridDegree %>% filter(Metric == "AllDegree") %>% 
  mutate(Degree = ifelse(Degree>300, 300, Degree)*Richness) %>%
  ggplot(aes(X, Y, fill = Degree, colour = Degree)) +
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log10(Richness+1))) +
  coord_fixed() + 
  guides(alpha = "none") +
  labs(fill = "All links", colour = "All links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[1]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[1]) +
  theme_void() + 
  theme(legend.position = "bottom",
        legend.text = element_text(hjust = 1, angle = 40))

Sum_Map_In <- GridDegree %>% filter(Metric == "InDegree")  %>% 
  mutate(Degree = Degree*Richness) %>%
  ggplot(aes(X, Y, fill = Degree, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log10(Richness+1))) +
  coord_fixed() + 
  guides(alpha = "none") +
  labs(fill = "Within-order links", colour = "Within-order links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2]) +
  theme_void() + 
  theme(legend.position = "bottom",
        legend.text = element_text(hjust = 1, angle = 40)) 

Sum_Map_Out <- GridDegree %>% filter(Metric == "OutDegree")  %>% 
  mutate(Degree = ifelse(Degree>150, 150, ifelse(Degree<0, 0, Degree))*Richness) %>%
  ggplot(aes(X, Y, fill = Degree, colour = Degree)) + 
  geom_tile(fill = "grey", colour = "grey") +
  geom_tile(aes(alpha = log10(Richness+1))) +
  coord_fixed() + 
  guides(alpha = "none") +
  labs(fill = "Out-of-order links", colour = "Out-of-order links") +
  scale_colour_continuous_sequential(palette = AlberPalettes[3]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[3]) +
  theme_void() + 
  theme(legend.position = "bottom",
        legend.text = element_text(hjust = 1, angle = 40)) 

TextSize = 3
AxisTextX = 8
AxisTextY = 10

Sum_Taxon_All <- 
  BarGraph(Panth1, "hOrder", "AllPredDegree", Just = T, Order = T, Text = "N", TextSize = TextSize, Fun = sum) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
  labs(x = NULL, y = "All links") +
  theme(axis.text.x = element_text(size = AxisTextX),
        axis.text.y = element_text(size = AxisTextY)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

Sum_Taxon_In <- BarGraph(Panth1, "hOrder", "InDegree", Just = T, Order = T, Text = "N", TextSize = TextSize, Fun = sum) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]]) +
  labs(x = NULL, y = "Within-order links") +
  theme(axis.text.x = element_text(size = AxisTextX),
        axis.text.y = element_text(size = AxisTextY)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

Sum_Taxon_Out <- BarGraph(Panth1, "hOrder", "OutDegree", Just = T, Order = T, Text = "N", TextSize = TextSize, Fun = sum) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[3]]) +
  labs(x = NULL, y = "Out-of-order links") +
  theme(axis.text.x = element_text(size = AxisTextX),
        axis.text.y = element_text(size = AxisTextY)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

Row_All <- plot_grid(Sum_Taxon_All, Sum_Map_All, 
                     nrow = 1, rel_widths = c(1.5,1), 
                     labels = LETTERS[1:2], axis = "t")

Row_In <- plot_grid(Sum_Taxon_In, Sum_Map_In, nrow = 1, rel_widths = c(1.5,1), 
                    labels = LETTERS[1:2+2], axis = "t")

Row_Out<- plot_grid(Sum_Taxon_Out, Sum_Map_Out, nrow = 1, rel_widths = c(1.5,1), 
                    labels = LETTERS[1:2+4], axis = "t")

SumPlot <- plot_grid(Row_All, Row_In, Row_Out, nrow = 3)

SumPlot %>% save_plot(filename = "SIFigures/Figure SI9.jpeg",
                      base_aspect_ratio = 1,
                      base_height = 9)

# Deviance contributions ####

Resps %>% lapply(function(a){ 
  
  gather(DevianceDFList[[a]], "Model", "Deviance", Model_Deviance, Total_Deviance) %>% 
    ggplot(aes(Model, Deviance, fill = Var)) + geom_col(position = "stack", colour = "black") + 
    lims(y = c(0,1)) +
    labs(fill = "Variable") +
    scale_fill_discrete_sequential(palette = "Plasma", rev = F,
                                   labels = c("Domestic", "Citations", "Space == 0", "Space", "Phylogeny", "Species")) +
    scale_x_discrete(labels = c("Model Deviance", "Total Deviance")) +
    ggtitle(a)
  
}) %>% plot_grid(plotlist = ., ncol = 5) %>% save_plot(file = "SIFigures/DevianceOutput.jpeg", 
                                                       base_aspect_ratio = 1, base_width = 15)


Resps %>% lapply(function(a){ 
  
  gather(DevianceDFList[[a]], "Model", "Deviance", Model_Deviance, Total_Deviance)
  
})  %>% bind_rows(.id = "Response") %>% mutate(Response = factor(RespLabels[as.numeric(Response)], levels = RespLabels)) %>%
  ggplot(aes(Model, Deviance, fill = Var)) + geom_col(position = "stack", colour = "black") + 
  lims(y = c(0,1)) +
  labs(fill = "Variable", x = NULL) +
  scale_fill_discrete_sequential(palette = "Plasma", rev = F,
                                 labels = c("Domestic", "Citations", "Space == 0", "Space", "Phylogeny", "Species")) +
  scale_x_discrete(labels = c("Model Deviance", "Total Deviance")) +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(hjust = 1, angle = 45)) + 
  facet_wrap(~Response, nrow = 1) +
  ggsave(file = "SIFigures/DevianceOutput.jpeg", units = "mm", height = 100, width = 200)

