# 0_Kludging spatial data import ####

# This was all done in 2 days and is probably circuitous
# Because I mashed everything into pixels and data frames rather than working
# with spatial dataframes :howdy: :grimace:

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

if(file.exists("data/MammalRanges.Rdata")) load(file = "data/MammalRanges.Rdata") else{
  
  mammal_shapes <- st_read("TERRESTRIAL_MAMMALS")
  
  #mammal_shapes <- st_transform(mammal_shapes, 54009) # Mollweide projection 
  mammal_shapes <- st_transform(mammal_shapes, 
                                "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
  
  # Mollweide projection = +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
  # This projection retains grid size as much as possible, but at the expense of shape
  
  mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
  mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
  mammal_shapes_red <- mammal_shapes[mammal_shapes$binomial%in%unique(Hosts$Sp),]
  mammal_raster <- raster(mammal_shapes_red, res = 50000) # NB units differ from Mercator!
  
  MammalRanges <- fasterize(mammal_shapes_red, mammal_raster, by = "binomial")
  #save(MammalRanges, file = "data/MammalRanges.Rdata")
  
  remove("mammal_raster", "mammal_shapes", "mammal_shapes_red")
  
}

LeftOut <- Hosts[!Hosts$Sp%in%names(MammalRanges),"Sp"]

if(file.exists("data/MammalRanges2.Rdata")) load(file = "data/MammalRanges2.Rdata") else{
  
  mammal_shapes2 <- st_read("Mammals_Terrestrial")
  mammal_shapes2 <- st_transform(mammal_shapes2, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection
  
  mammal_shapes2$binomial = str_replace(mammal_shapes2$BINOMIAL, " ", "_")
  mammal_shapes2 <- mammal_shapes2[order(mammal_shapes2$binomial),]
  mammal_shapes_red2 <- mammal_shapes2[mammal_shapes2$binomial%in%LeftOut,]
  
  mammal_raster2 <- raster(mammal_shapes_red2, res = 50000) # NB units differ from Mercator!
  
  MammalRanges2 <- fasterize(mammal_shapes_red2, mammal_raster2, by = "binomial")
  
  remove("mammal_raster2", "mammal_shapes2", "mammal_shapes_red2")
  
}

data("wrld_simpl")
WorldMap <- spTransform(wrld_simpl, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# I HATE this 

WorldPolygons <- lapply(WorldMap@polygons, 
                        function(b){ lapply(b@Polygons, 
                                            function(a) as.data.frame(a@coords))}) %>% 
  unlist(recursive = F)


WorldPolygons <- lapply(1:length(WorldMap@polygons), 
                        function(b){ lapply(WorldMap@polygons[[b]]@Polygons, 
                                            function(a){ as.data.frame(a@coords) %>%
                                                mutate(Country = WorldMap@data[b,"NAME"])})}) %>% 
  unlist(recursive = F)

for(x in 1:length(WorldPolygons)) WorldPolygons[[x]]$group <- x 

WorldMap <- bind_rows(WorldPolygons) %>% dplyr::rename(long = V1, lat = V2)

WorldMap <- WorldMap %>% group_by(Country, group) %>%  # getting rid of small places 
  summarise(N = n()) %>% filter(N>10) %>% inner_join(WorldMap)

CountryCentroids <- WorldMap %>% group_by(Country) %>% 
  summarise(long = mean(long), lat = mean(lat)) %>% data.frame

save(WorldMap, file = "Output Files/WorldMap.Rdata")

# THIS DATA FRAME TAKES A LOT OF MEMORY - convert to sparse matrix 
# Or learn more raster methods before pub ####

if(!file.exists("data/RangeOverlap.Rdata")){
  
  Valuedf <- data.frame(getValues(MammalRanges))
  Valuedf2 <- reshape2::melt(Valuedf)
  Valuedf2$x <- rep(1:MammalRanges[[1]]@ncols, MammalRanges[[1]]@nrows)
  Valuedf2$y <- rep(MammalRanges[[1]]@nrows:1, each = MammalRanges[[1]]@ncols)
  Valuedf2 <- Valuedf2 %>% slice(-which(is.na(value)))
  
  Valuedf3 <- data.frame(getValues(MammalRanges2))
  Valuedf4 <- reshape2::melt(Valuedf3)
  Valuedf4$x <- rep(1:MammalRanges2[[1]]@ncols, MammalRanges2[[1]]@nrows)
  Valuedf4$y <- rep(MammalRanges2[[1]]@nrows:1, each = MammalRanges2[[1]]@ncols)
  Valuedf4 <- Valuedf4 %>% slice(-which(is.na(value)))
  
  Rangedf <- rbind(Valuedf2,Valuedf4) %>% 
    dplyr::rename(Host = variable, Presence = value)
  
  Rangedf$GridID <- with(Rangedf, paste(x, y))
  
  Range0 <- levels(Rangedf$Host)[which(table(Rangedf$Host)==0)] # Hosts that have no spatial records??
  Rangedf <- droplevels(Rangedf) 
  Rangedf <- Rangedf[order(Rangedf$Host),]
  
  # Using igraph to project it
  
  RangeOverlap <- matrix(0, nrow = nlevels(Rangedf$Host), ncol = nlevels(Rangedf$Host))
  dimnames(RangeOverlap) <- list(levels(Rangedf$Host),levels(Rangedf$Host))
  
  for(x in levels(Rangedf$Host)){
    
    if(x == first(levels(Rangedf$Host))) t1 <- Sys.time()
    
    Grids <- Rangedf[Rangedf$Host==x,"GridID"]
    SubRangedf <- Rangedf[Rangedf$GridID %in% Grids,]
    
    RangeOverlap[x,] <- table(SubRangedf$Host)
    
    print(x)
    
    if(x == last(levels(Rangedf$Host))) t2 <- Sys.time()
    
  }
  
  RangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)), nrow(RangeOverlap))
  RangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)), nrow(RangeOverlap))
  
  RangeAdj <- RangeOverlap/(RangeA + RangeB - RangeOverlap) # Weighted evenly
  
  save(RangeAdj, file = "data/RangeOverlap.Rdata")
  
  remove(MammalRanges, MammalRanges2)
  
} else load("data/RangeOverlap.Rdata")

Hosts[Hosts$Sp%in%rownames(RangeAdj),"S.Greg1"] <- rowSums(RangeAdj[as.character(Hosts[Hosts$Sp%in%rownames(RangeAdj),"Sp"]),])

# Making polygons for display ####

if(file.exists("data/HostPolygons.Rdata")) load("data/HostPolygons.Rdata") else{
  
  HostPolygons <- lapply(levels(Valuedf2$variable), function(x) {
    
    if(!x%in%Range0){
      
      r <- MammalRanges[[x]] > -Inf
      
      r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
        mutate(Host = x) %>% return
    }
    
  }) %>% bind_rows()
  
  HostPolygons2 <- lapply(levels(Valuedf4$variable), function(x) {
    
    if(!x%in%Range0){
      
      r <- MammalRanges2[[x]] > -Inf
      
      r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
        mutate(Host = x) %>% return
    }
    
  }) %>% bind_rows()
  
  HostPolygons <- bind_rows(HostPolygons, HostPolygons2)
  
  save(HostPolygons, file = "data/HostPolygons.Rdata")}

# Deriving Host Centroids and Merging ####

HostCentroids <- data.frame(LongMean = with(HostPolygons, tapply(long, Host, mean)),
                            LatMean = with(HostPolygons, tapply(lat, Host, mean)),
                            Host = unique(HostPolygons$Host))

Hosts <- merge(Hosts, HostCentroids, all.x = T, by.x = "Sp", by.y = "Host")

# Viral spatial data adapted from spatial data for their hosts ####

# Making Viral Associations and Polygons ####

VirusAssocs <- apply(M, 1, function(a) names(a[a>0]))

load("~/Albersnet/data/FullPolygons.Rdata")

CodeRoot <- "R Code/0_Data Import"

if(file.exists("~/Albersnet/data/FullRangeOverlap.Rdata")) load("~/Albersnet/data/FullRangeOverlap.Rdata") else{
  source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
}

detach(package:raster)
