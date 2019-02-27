
# Summarising predicted degree by grid 

library(tidyverse); library(raster); library(colorspace)

load("~/Albersnet/data/FullMammalRanges.Rdata")
load("Output Files/Panth1.Rdata")

FullValuedf <- data.frame(getValues(FullMammalRanges)); detach(package:raster)
FullValuedf2 <- reshape2::melt(FullValuedf)
FullValuedf2$x <- rep(1:FullMammalRanges[[1]]@ncols, FullMammalRanges[[1]]@nrows)
FullValuedf2$y <- rep(FullMammalRanges[[1]]@nrows:1, each = FullMammalRanges[[1]]@ncols)

detach(package:raster)

FullRangedf <- FullValuedf2 %>% 
  filter(!is.na(value)) %>% droplevels %>%
  dplyr::rename(Host = variable, Presence = value)

save(FullRangedf, file = "Output Files/FullRangedf.Rdata")
#load("Output Files/FullRangedf.Rdata")

FullRangedf2 <- FullRangedf %>% left_join(Panth1[,c("Sp","AllPredDegree", "InDegree", "OutDegree")], 
                                          by = c("Host" = "Sp")) %>%
  filter(!Host=="Ursus_maritimus")

GridDegree <- FullRangedf2 %>% group_by(x,y) %>% 
  summarise_at(vars(ends_with("Degree")), function(a) mean(a, na.rm = T))

#save(GridDegree, file = "Output Files/GridDegree.Rdata")

GridDegree2 <- gather(GridDegree, key = "Metric", value = "Degree", ends_with("Degree"))

save(GridDegree2, file = "Output Files/GridDegree2.Rdata")

GridDegree2 %>% filter(Metric == "AllPredDegree") %>% filter(abs(scale(Degree))<310) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(AllPredDegree = "All Links"))) +
  coord_fixed() +  
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2])

GridDegree2 %>% filter(Metric == "InDegree") %>% filter(Degree<180) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(InDegree = "Within-Order Links"))) +
  coord_fixed() +  
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[2]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[2])

GridDegree2 %>% filter(Metric == "OutDegree") %>% filter(Degree<190) %>%
  ggplot(aes(x, y, fill = Degree, colour = Degree)) + geom_tile() +
  facet_wrap(~Metric, nrow = 3, labeller = labeller(Metric = c(OutDegree = "Out-of-Order Links"))) +
  coord_fixed() +  
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous_sequential(palette = AlberPalettes[3]) +  
  scale_fill_continuous_sequential(palette = AlberPalettes[3])
