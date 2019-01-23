
# Data check: what mammals don't have spatial data? 

LeftOut <- Hosts[is.na(Hosts$LongMean),"Sp"] %>% as.character

mammal_shapes2 <- st_read("Mammals_Terrestrial")
mammal_shapes2 <- st_transform(mammal_shapes2, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection

mammal_shapes2$binomial = str_replace(mammal_shapes2$BINOMIAL, " ", "_")
mammal_shapes2 <- mammal_shapes2[order(mammal_shapes2$binomial),]
mammal_shapes_red2 <- mammal_shapes2[mammal_shapes2$binomial%in%LeftOut,]

mammal_raster2 <- raster(mammal_shapes_red2, res = 50000) # NB units differ from Mercator!

MammalRanges2 <- fasterize(mammal_shapes_red2, mammal_raster2, by = "binomial")

LeftOut%in%names(MammalRanges2)

NeverFinding <- setdiff(LeftOut, c(mammal_shapes$binomial)) # All domestic or marine ####
ToRemedy <- intersect(LeftOut, c(mammal_shapes$binomial))

Valuedf <- data.frame(getValues(MammalRanges))

#SAVEHostsLeftout <- levels(Rangedf$Host)[which(table(Rangedf$Host)==0)]

levels(Rangedf$Host)[which(table(Rangedf$Host)==0)]
