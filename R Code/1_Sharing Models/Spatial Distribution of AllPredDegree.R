
# Summarising predicted degree by grid 

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools)

load("~/Albersnet/data/FullMammalRanges.Rdata")
load("~/Albersnet/data/FullMammalRanges2.Rdata")

<<<<<<< HEAD

=======
<<<<<<< HEAD
=======
>>>>>>> a4cef8a39db290aea40687955d8fdd8fcded2a7e
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

with(FullRangedf, )


<<<<<<< HEAD
=======
>>>>>>> f2e086d8c659c0f652c2f24c7a47811ae8c4372e
>>>>>>> a4cef8a39db290aea40687955d8fdd8fcded2a7e
# Making shapefile ####

# 1. Put degree, etc. into shapefile attribute table


# 2. Convert spatial polygons data frame into sf object

ms <- st_as_sf(
  FullPolygons, 
  coords = c('long', 'lat')
)

nc_sp <- sf:::as_Spatial(FullPolygons)

maptools::unionSpatialPolygons(SpP, IDs)

FullPolygons %>% group_by(Host, group) #%>% unlist()

# trying the sf I already have ####

<<<<<<< HEAD
load("~/Albersnet/AllPredDegree.Rdata")

mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
mammal_shapes$AllPredDegree <- AllPredDegrees[mammal_shapes$binomial]


mammal_shapes <- st_read("~/Albersnet shapefiles/TERRESTRIAL_MAMMALS (new)")


=======
<<<<<<< HEAD
load("~/Albersnet/AllPredDegree.Rdata")
=======
mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
mammal_shapes$AllPredDegree <- AllPredDegrees[mammal_shapes$binomial]
>>>>>>> f2e086d8c659c0f652c2f24c7a47811ae8c4372e

mammal_shapes <- st_read("~/Albersnet shapefiles/TERRESTRIAL_MAMMALS (new)")

<<<<<<< HEAD
>>>>>>> a4cef8a39db290aea40687955d8fdd8fcded2a7e
#mammal_shapes <- st_transform(mammal_shapes, 54009) # Mollweide projection 
mammal_shapes <- st_transform(mammal_shapes, 
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 

# Mollweide projection = +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
# This projection retains grid size as much as possible, but at the expense of shape

mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
mammal_shapes$AllPredDegree <- AllPredDegrees[mammal_shapes$binomial]

mammal_raster_full <- raster(mammal_shapes, res = 50000) # NB units differ from Mercator!

<<<<<<< HEAD

=======
=======
>>>>>>> f2e086d8c659c0f652c2f24c7a47811ae8c4372e
>>>>>>> a4cef8a39db290aea40687955d8fdd8fcded2a7e
DegreeRanges <- fasterize(mammal_shapes,
                          mammal_raster_full, 
                          by = "binomial", 
                          field = "AllPredDegree")

DegreeRanges2 <- raster::stackApply(DegreeRanges, 3646, fun = "mean", na.rm = T)
<<<<<<< HEAD


DegreeRanges2 <- raster::calc(DegreeRanges, fun = "mean")
DegreeRanges2 <- raster::calc(DegreeRanges@data@values, fun = "mean")

=======
<<<<<<< HEAD

DegreeRanges2 <- raster::calc(DegreeRanges, fun = "mean")
DegreeRanges2 <- raster::calc(DegreeRanges@data@values, fun = "mean")
=======
>>>>>>> f2e086d8c659c0f652c2f24c7a47811ae8c4372e
>>>>>>> a4cef8a39db290aea40687955d8fdd8fcded2a7e

DegreeRanges2 <- raster::calc(DegreeRanges, fun = "mean")
DegreeRanges2 <- raster::calc(DegreeRanges@data@values, fun = "mean")


# 3. Use Noam's fasterize package to collapse polygons into a raster
# with function='mean' acting on whatever variable you care about




