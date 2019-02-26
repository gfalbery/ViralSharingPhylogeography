
# Summarising predicted degree by grid 

library(tidyverse); library(raster)

load("~/Albersnet/data/FullMammalRanges.Rdata")

FullValuedf <- data.frame(getValues(FullMammalRanges))
FullValuedf2 <- reshape2::melt(FullValuedf)
FullValuedf2$x <- rep(1:FullMammalRanges[[1]]@ncols, FullMammalRanges[[1]]@nrows)
FullValuedf2$y <- rep(FullMammalRanges[[1]]@nrows:1, each = FullMammalRanges[[1]]@ncols)

FullRangedf <- FullValuedf2 %>% 
  dplyr::rename(Host = variable, Presence = value)

FullRangedf$GridID <- with(FullRangedf, paste(x, y))

FullRangedf <- droplevels(FullRangedf) 
FullRangedf <- FullRangedf[order(FullRangedf$Host),]

FullRangedf <- FullRangedf %>% left_join(Panth1, by.x = "Host", by.y = "Sp", all.x = T)

GridDegree <- with(FullRangedf, tapply(AllPredDegree, list(x, y), function(a) mean(a, na.rm = T))) %>% 
  reshape2::melt() %>% na.omit %>% dplyr::rename(AllDegree = value)

GridDegree$InDegree <- with(FullRangedf, tapply(InDegree, list(x, y), function(a) mean(a, na.rm = T))) %>% 
  reshape2::melt() %>% na.omit %>% select(value)

GridDegree$OutDegree <- with(FullRangedf, tapply(OutDegree, list(x, y), function(a) mean(a, na.rm = T))) %>% 
  reshape2::melt() %>% na.omit %>% select(value)

GridDegree2 <- gather(GridDegree, key = "Metric", value = "Degree", AllDegree, InDegree, OutDegree)

ggplot(GridDegree2, aes(x, y, fill = Degree)) +
  facet_wrap(~Metric) +
  geom_tile() + coord_fixed() + 
  scale_fill_gradientn(colours = c(AlberColours[1], 
                                   AlberColours[1], "white",
                                   AlberColours[2], AlberColours[2]), values = c(0,0.33,0.8,1))

save(GridDegree, file = "Output Files/GridDegree.Rdata")
save(GridDegree2, file = "Output Files/GridDegree2.Rdata")

