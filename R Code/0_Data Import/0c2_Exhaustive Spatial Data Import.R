# Creating exhaustive mammal spatial data ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools)

# Importing/making ranges ####

if(file.exists("data/FullMammalRanges.Rdata")) load("data/FullMammalRanges.Rdata") else{
  
  mammal_shapes <- st_read("~/Albersnet shapefiles/TERRESTRIAL_MAMMALS (new)")
  
  #mammal_shapes <- st_transform(mammal_shapes, 54009) # Mollweide projection 
  mammal_shapes <- st_transform(mammal_shapes, 
                                "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
  
  # Mollweide projection = +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
  # This projection retains grid size as much as possible, but at the expense of shape
  
  mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
  mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
  mammal_raster_full <- raster(mammal_shapes, res = 50000) # NB units differ from Mercator!
  
  FullMammalRanges <- fasterize(mammal_shapes, mammal_raster_full, by = "binomial")
  save(FullMammalRanges, file = "data/FullMammalRanges.Rdata")
  
}

# Trying earlier dataset ####

if(file.exists("data/FullMammalRanges2.Rdata")) load("data/FullMammalRanges2.Rdata") else{
  
  mammal_shapes2 <- st_read("~/Albersnet shapefiles/Mammals_Terrestrial (old)")
  mammal_shapes2 <- st_transform(mammal_shapes2, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection
  
  mammal_shapes2$binomial = str_replace(mammal_shapes2$BINOMIAL, " ", "_")
  mammal_shapes2 <- mammal_shapes2[order(mammal_shapes2$binomial),]
  mammal_shapes_red2 <- mammal_shapes2[!mammal_shapes2$binomial%in%names(FullMammalRanges),]
  
  #mammal_raster_full2 <- raster(mammal_shapes2, res = 50000) # NB units differ from Mercator!
  
  FullMammalRanges2 <- fasterize(mammal_shapes_red2, mammal_raster_full2, by = "binomial")
  
  save(FullMammalRanges2, file = "data/FullMammalRanges2.Rdata")
  
}

# Converting these to meaningful values ####

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

# Could use igraph to project it into bipartite host-grid matrix
# Or could do this bullshit

FullRangeOverlap <- matrix(0, nrow = nlevels(FullRangedf$Host), ncol = nlevels(FullRangedf$Host))
dimnames(FullRangeOverlap) <- list(levels(FullRangedf$Host),levels(FullRangedf$Host))

for(x in levels(FullRangedf$Host)){
  
  if(x == first(levels(FullRangedf$Host))) t1 <- Sys.time()
  
  Grids <- FullRangedf[FullRangedf$Host==x,"GridID"]
  SubFullRangedf <- FullRangedf[FullRangedf$GridID %in% Grids,]
  
  FullRangeOverlap[x,] <- table(SubFullRangedf$Host)
  
  print(x)
  
  if(x == last(levels(FullRangedf$Host))) t2 <- Sys.time()
  
}

FullRangeA = matrix(rep(diag(FullRangeOverlap), nrow(FullRangeOverlap)), nrow(FullRangeOverlap))
FullRangeB = matrix(rep(diag(FullRangeOverlap), each = nrow(FullRangeOverlap)), nrow(FullRangeOverlap))

FullRangeAdj1 <- FullRangeOverlap/(FullRangeA + FullRangeB - FullRangeOverlap) # Weighted evenly
FullRangeAdj2 <- FullRangeOverlap/(FullRangeA) # Asymmetrical

save(FullRangeAdj1, file = "data/FullRangeOverlap.Rdata")

# Making polygons for display ####

FullPolygons <- lapply(levels(FullValuedf2$variable), function(x) {
  
  if(!x%in%Range0){
    
    r <- FullMammalRanges[[x]] > -Inf
    
    r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
      mutate(Host = x) %>% return
  }
  
}) %>% bind_rows()

FullPolygons2 <- lapply(levels(FullValuedf4$variable), function(x) {
  
  if(!x%in%Range0){
    
    r <- FullMammalRanges2[[x]] > -Inf
    
    r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
      mutate(Host = x) %>% return
  }
  
}) %>% bind_rows()

FullPolygons <- bind_rows(FullPolygons, FullPolygons2)

save(FullPolygons, file = "data/FullPolygons.Rdata")
