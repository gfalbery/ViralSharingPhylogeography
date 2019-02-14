
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

FullRangedf$AllPredDegree <- AllPredDegrees[as.character(FullRangedf$Host)]

head(FullRangedf)
