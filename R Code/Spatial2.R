# Rscript "R Code/Spatial2.R"

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

#load("~/LargeFiles/OverList1.Rdata")
#print(dim(OverList1))

#load("~/LargeFiles/OverList2.Rdata")
#print(dim(OverList2))

if(file.exists("~/LargeFiles/FullRangedf.Rdata")){load("~/LargeFiles/FullRangedf.Rdata"); print("Loaded!") }else{
  
  FullRangedf <- rbind(OverList1 %>% bind_rows() %>% mutate(Pass = 1), OverList2 %>% bind_rows() %>% mutate(Pass = 2))
  
  FullRangedf <- FullRangedf %>% slice(order(Pass, Species))
  
  yDiff <- (FullRangedf %>% group_by(Pass) %>% summarise(Diff = max(y)))$Diff %>% diff
  
  FullRangedf <- FullRangedf %>% mutate(y = ifelse(Pass==1, y + yDiff, y))
  
  Pass1Sp <- FullRangedf %>% filter(Pass == 1)
  Pass2Sp <- FullRangedf %>% filter(Pass == 2)
  
  SecondSp <- setdiff(Pass2Sp$Species, Pass1Sp$Species)
  
  FullRangedf <- FullRangedf %>% filter(Pass==1|Species %in% SecondSp)
  
  save(FullRangedf, file = "~/LargeFiles/FullRangedf.Rdata")
  
  print("Saved!")
}

load("Output Files/Panth1.Rdata")

print("Joining rangedf and Panth1")

FullRangedf2 <- FullRangedf %>% #filter(Pass==1|Species %in% SecondSp) %>%
  left_join(Panth1[,c("Sp","AllPredDegree", "InDegree", "OutDegree")], 
                                          by = c("Species" = "Sp")) %>%
  filter(!Species=="Ursus_maritimus")

GridDegree <- FullRangedf2 %>% group_by(x,y) %>% 
  summarise_at(vars(ends_with("Degree")), function(a) mean(a, na.rm = T))

PopDegree <- FullRangedf2 %>% group_by(x,y) %>% 
  summarise(Density = n())

#save(GridDegree, file = "Output Files/GridDegree.Rdata")

GridDegree2 <- gather(GridDegree, key = "Metric", value = "Degree", ends_with("Degree")) %>% 
  left_join(PopDegree, by = c("x","y"))

save(GridDegree2, file = "~/LargeFiles/GridDegree2.Rdata")

