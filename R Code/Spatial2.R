# Rscript "R Code/Spatial2.R"

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

load("~/LargeFiles/OverList1.Rdata")
print(dim(OverList1))

load("~/LargeFiles/OverList2.Rdata")
print(dim(OverList2))

FullRangedf <- rbind(OverList1 %>% bind_rows() %>% mutate(Pass = 1), OverList2 %>% bind_rows() %>% mutate(Pass = 2))

FullRangedf <- FullRangedf %>% slice(order(Pass, Host))

save(FullRangedf, file = "~/LargeFiles/FullRangedf.Rdata")

print("Saved!")

load("Output Files/Panth1.Rdata")

print("Joining rangedf and Panth1")

FullRangedf2 <- FullRangedf %>% left_join(Panth1[,c("Sp","AllPredDegree", "InDegree", "OutDegree")], 
                                          by = c("Host" = "Sp")) %>%
  filter(!Host=="Ursus_maritimus")

GridDegree <- FullRangedf2 %>% group_by(x,y) %>% 
  summarise_at(vars(ends_with("Degree")), function(a) mean(a, na.rm = T))

PopDegree <- FullRangedf2 %>% group_by(x,y) %>% 
  summarise(Density = n())

#save(GridDegree, file = "Output Files/GridDegree.Rdata")

GridDegree2 <- gather(GridDegree, key = "Metric", value = "Degree", ends_with("Degree")) %>% 
  left_join(PopDegree, by = c("x","y"))

save(GridDegree2, file = "Large Files/GridDegree2.Rdata")

