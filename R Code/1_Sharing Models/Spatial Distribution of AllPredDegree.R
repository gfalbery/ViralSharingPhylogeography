
# Summarising predicted degree by grid 

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools)

load("~/Albersnet/data/FullMammalRanges.Rdata")
load("~/Albersnet/data/FullMammalRanges2.Rdata")

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

load("~/Albersnet/AllPredDegree.Rdata")

mammal_shapes <- st_read("~/Albersnet shapefiles/TERRESTRIAL_MAMMALS (new)")

#mammal_shapes <- st_transform(mammal_shapes, 54009) # Mollweide projection 
mammal_shapes <- st_transform(mammal_shapes, 
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 

# Mollweide projection = +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
# This projection retains grid size as much as possible, but at the expense of shape

mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
mammal_shapes$AllPredDegree <- AllPredDegrees[mammal_shapes$binomial]

mammal_raster_full <- raster(mammal_shapes, res = 50000) # NB units differ from Mercator!

DegreeRanges <- fasterize(mammal_shapes,
                          mammal_raster_full, 
                          by = "binomial", 
                          field = "AllPredDegree")

DegreeRanges2 <- raster::stackApply(DegreeRanges, 3646, fun = "mean", na.rm = T)

DegreeRanges2 <- raster::calc(DegreeRanges, fun = "mean")
DegreeRanges2 <- raster::calc(DegreeRanges@data@values, fun = "mean")

plot(DegreeRanges2)


# 3. Use Noam's fasterize package to collapse polygons into a raster
# with function='mean' acting on whatever variable you care about




