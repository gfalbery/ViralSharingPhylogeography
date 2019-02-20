
# Summarising predicted degree by grid 

library(tidyverse)

load("~/Albersnet/data/FullMammalRanges.Rdata")
load("~/Albersnet/data/FullMammalRanges2.Rdata")

FullValuedf <- data.frame(getValues(FullMammalRanges))
FullValuedf2 <- reshape2::melt(FullValuedf)
FullValuedf2$x <- rep(1:FullMammalRanges[[1]]@ncols, FullMammalRanges[[1]]@nrows)
FullValuedf2$y <- rep(FullMammalRanges[[1]]@nrows:1, each = FullMammalRanges[[1]]@ncols)

FullValuedf3 <- data.frame(getValues(FullMammalRanges2))
FullValuedf4 <- reshape2::melt(FullValuedf3)
FullValuedf4$x <- rep(1:FullMammalRanges2[[1]]@ncols, FullMammalRanges2[[1]]@nrows)
FullValuedf4$y <- rep(FullMammalRanges2[[1]]@nrows:1, each = FullMammalRanges2[[1]]@ncols)

FullRangedf <- rbind(FullValuedf2[!is.na(FullValuedf2$value),],FullValuedf4[!is.na(FullValuedf4$value),]) # This is where a load of them were lost ####
FullRangedf <- FullRangedf %>% 
  dplyr::rename(Host = variable, Presence = value)

FullRangedf$GridID <- with(FullRangedf, paste(x, y))

Range0 <- levels(FullRangedf$Host)[which(table(FullRangedf$Host)==0)] # Hosts that have no spatial records??
FullRangedf <- droplevels(FullRangedf) 
FullRangedf <- FullRangedf[order(FullRangedf$Host),]

Incorporated <- c(names(FullMammalRanges), names(FullMammalRanges2))

FullRangedf$AllPredDegree <- NA

for(x in 0:400){
  if(x %% 100==0) print(x)
  FullRangedf$AllPredDegree[1:10000+(10000*x)] <- AllPredDegrees[as.character(FullRangedf$Host)[1:10000+(10000*x)]]
}

FullRangedf$AllPredDegree[(10000*(x-1)):nrow(FullRangedf)] <- AllPredDegrees[as.character(FullRangedf$Host)[1+(10000*(x-1)):nrow(FullRangedf)]]

AllDF <- data.frame(Host = names(AllPredDegrees), AllPredDegree = AllPredDegrees)

FullRangedf <- FullRangedf %>% left_join()

OutDegrees <- Panth1$OutDegree
names(OutDegrees) <- Panth1$Sp

FullRangedf$OutDegree <- NA

for(x in 0:400){
  if(x %% 100==0) print(x)
  FullRangedf$OutDegree[1:10000+(10000*x)] <- OutDegrees[as.character(FullRangedf$Host)[1:10000+(10000*x)]]
}

FullRangedf$OutDegree[(10000*(x-1)):nrow(FullRangedf)] <- OutDegrees[as.character(FullRangedf$Host)[(10000*(x-1)):nrow(FullRangedf)]]

InDegrees <- Panth1$InDegree
names(InDegrees) <- Panth1$Sp

FullRangedf$InDegree <- NA

for(x in 0:400){
  if(x %% 100==0) print(x)
  FullRangedf$InDegree[1:10000+(10000*x)] <- InDegrees[as.character(FullRangedf$Host)[1:10000+(10000*x)]]
}

FullRangedf$InDegree[(10000*(x-1)):nrow(FullRangedf)] <- InDegrees[as.character(FullRangedf$Host)[(10000*(x-1)):nrow(FullRangedf)]]

unique(FullRangedf$Host)

FullRangedf$Pass <- ifelse(FullRangedf$Host%in%FullValuedf2$variable,1,2)

GridDegree <- with(FullRangedf[FullRangedf$Pass==1,], tapply(AllPredDegree, list(x, y), function(a) mean(a, na.rm = T))) %>% 
  reshape2::melt() %>% na.omit %>% dplyr::rename(AllDegree = value)

GridDegree$InDegree <- with(FullRangedf[FullRangedf$Pass==1,], tapply(InDegree, list(x, y), function(a) mean(a, na.rm = T))) %>% 
  reshape2::melt() %>% na.omit %>% select(value)

GridDegree$OutDegree <- with(FullRangedf[FullRangedf$Pass==1,], tapply(OutDegree, list(x, y), function(a) mean(a, na.rm = T))) %>% 
  reshape2::melt() %>% na.omit %>% select(value)

GridDegree2 <- gather(GridDegree, key = "Metric", value = "Degree", AllDegree, InDegree, OutDegree)

ggplot(GridDegree, aes(Var1, Var2, fill = value)) + 
  geom_tile() + coord_fixed() + 
  scale_fill_gradientn(colours = c(AlberColours[1], 
                                   AlberColours[1], "white",
                                   AlberColours[2], AlberColours[2]), values = c(0,0.33,0.8,1))
