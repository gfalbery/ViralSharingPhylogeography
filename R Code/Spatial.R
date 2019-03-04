# Creating exhaustive mammal spatial data ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

# Importing/making ranges ####]

print("Mammal Ranges 1!")

if(file.exists("~/LargeFiles/MammalRanges1.Rdata")) load("~/LargeFiles/MammalRanges1.Rdata") else{
  
  mammal_shapes <- st_read("~/ShapeFiles")
  
  mammal_shapes <- st_transform(mammal_shapes, 
                                "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
  
  mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
  mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
  mammal_raster_full <- raster(mammal_shapes, res = 20000)
  
  print("Fasterising!")
  
  FullMammalRanges <- fasterize(mammal_shapes, mammal_raster_full, by = "binomial")
  save(FullMammalRanges, file = "~/LargeFiles/MammalRanges1.Rdata")
  
}

# Trying earlier dataset ####

if(file.exists("~/LargeFiles/MammalRanges2.Rdata")) load("~/LargeFiles/MammalRanges2.Rdata") else{
  
  print("Mammal Ranges 2!")
  
  mammal_shapes2 <- st_read("~/ShapeFiles2")
  mammal_shapes2 <- st_transform(mammal_shapes2, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection
  
  mammal_shapes2$binomial = str_replace(mammal_shapes2$BINOMIAL, " ", "_")
  mammal_shapes2 <- mammal_shapes2[order(mammal_shapes2$binomial),]
  #mammal_shapes_red2 <- mammal_shapes2[!mammal_shapes2$binomial%in%names(FullMammalRanges),]
  
  mammal_raster_full2 <- raster(mammal_shapes2, res = 20000) # NB units differ from Mercator!
  
  print("Fasterising!")
  
  FullMammalRanges2 <- fasterize(mammal_shapes2, mammal_raster_full2, by = "binomial")
  
  save(FullMammalRanges2, file = "~/LargeFiles/MammalRanges2.Rdata")
  
}

# Converting these to meaningful values ####

print("Converting to values!")

if(file.exists("~/LargeFiles/FullValuedf1.Rdata")) load("~/LargeFiles/FullValuedf1.Rdata") else{
  
  FullValuedf <- data.frame(getValues(FullMammalRanges))
  
  print("Saving!")
  
  save(FullValuedf, file = "FullValuedf1.Rdata")
}

if(file.exists("~/LargeFiles/OverList1.Rdata")) load("~/LargeFiles/OverList1.Rdata") else{
  
  print("Making Overlist 1!")
  
  OverList1 <- lapply(1:ncol(FullValuedf), function(a){
    
    print(names(FullValuedf)[a])
    
    data.frame(Species = as.character(names(FullValuedf)[a]),
               Presence = FullValuedf[,a],
               x = rep(1:FullMammalRanges[[1]]@ncols, FullMammalRanges[[1]]@nrows),
               y = rep(FullMammalRanges[[1]]@nrows:1, each = FullMammalRanges[[1]]@ncols)
               
    ) %>% na.omit() %>% select(-Presence)
    
  }) #%>% bind_rows %>% mutate(Pass = 1)
  
  save(OverList1, file = "~/LargeFiles/OverList1.Rdata")
}

#remove(FullMammalRanges)
#remove(FullValuedf)

#FullValuedf2 <- reshape2::melt(FullValuedf)
#FullValuedf2$x <- rep(1:FullMammalRanges[[1]]@ncols, FullMammalRanges[[1]]@nrows)
#FullValuedf2$y <- rep(FullMammalRanges[[1]]@nrows:1, each = FullMammalRanges[[1]]@ncols)
#FullValuedf2$Pass <- 1

print("Converting the second one to values!")

if(file.exists("~/LargeFiles/FullValuedf3.Rdata")) load("~/LargeFiles/FullValuedf3.Rdata") else{
  
  FullValuedf3 <- data.frame(getValues(FullMammalRanges2))
  
  print("Saving!")
  
  save(FullValuedf3, file = "~/LargeFiles/FullValuedf3.Rdata")
}

if(file.exists("~/LargeFiles/OverList2.Rdata")) load("~/LargeFiles/OverList2.Rdata") else{
  
  print("Making Overlist 2!")
  
  OverList2 <- lapply(1:ncol(FullValuedf3), function(a){
    
    print(names(FullValuedf3)[a])
    
    data.frame(Species = as.character(names(FullValuedf3)[a]),
               Presence = FullValuedf3[,a],
               x = rep(1:FullMammalRanges2[[1]]@ncols, FullMammalRanges2[[1]]@nrows),
               y = rep(FullMammalRanges2[[1]]@nrows:1, each = FullMammalRanges2[[1]]@ncols)
               
    ) %>% na.omit() %>% select(-Presence)
    
  }) #%>% bind_rows %>% mutate(Pass = 2)
  
  save(OverList2, file = "~/LargeFiles/OverList2.Rdata")
  
}

print("Doing Overlap!")

if(file.exists("~/LargeFiles/FullRangeOverlap.Rdata")) load("~/LargeFiles/FullRangeOverlap.Rdata") else{
  
  FullRangeAdj1 <- PairsWisely(FullMammalRanges)
  
  save(FullRangeAdj1, file = "~/LargeFiles/FullRangeOverlap.Rdata")
  
}

#FullValuedf4 <- reshape2::melt(FullValuedf3)
#FullValuedf4$x <- rep(1:FullMammalRanges2[[1]]@ncols, FullMammalRanges2[[1]]@nrows)
#FullValuedf4$y <- rep(FullMammalRanges2[[1]]@nrows:1, each = FullMammalRanges2[[1]]@ncols)
#FullValuedf4$Pass <- 0

FullRangedf <- rbind(OverList1, OverList2) %>% # This is where a load of them were lost #### %>% 
  filter(!is.na(value)) %>% droplevels %>%
  dplyr::rename(Host = variable, Presence = value)

FullRangedf$GridID <- with(FullRangedf, paste(x, y))

FullRangedf <- FullRangedf[order(FullRangedf$Pass, FullRangedf$Host),]

save(FullRangedf, file = "~/LargeFiles/FullRangedf.Rdata")


# Could use igraph to project it into bipartite host-grid matrix
# Or could do this bullshit

# Making polygons for display ####

FullPolygons <- lapply(levels(FullValuedf2$variable), function(x) {
  
  if(!x%in%Range0){
    
    r <- FullMammalRanges[[x]] > -Inf
    
    r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
      mutate(Host = x) %>% return
  }
  
}) %>% bind_rows()

#FullPolygons2 <- lapply(levels(FullValuedf4$variable), function(x) {
#  
#  if(!x%in%Range0){
#    
#    r <- FullMammalRanges2[[x]] > -Inf
#    
#    r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
#      mutate(Host = x) %>% return
#  }
#  
#}) %>% bind_rows()

FullPolygons <- bind_rows(FullPolygons, FullPolygons2)

save(FullPolygons, file = "~/LargeFiles/FullPolygons.Rdata")

detach(package:raster)

