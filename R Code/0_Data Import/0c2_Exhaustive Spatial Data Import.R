# Creating exhaustive mammal spatial data ####

# Rscript "R Code/Spatial.R"

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
  mammal_raster_full <- raster(mammal_shapes, res = 25000)
  
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
  
  mammal_raster_full2 <- raster(mammal_shapes2, res = 25000) # NB units differ from Mercator!
  
  print("Fasterising!")
  
  FullMammalRanges2 <- fasterize(mammal_shapes2, mammal_raster_full2, by = "binomial")
  
  save(FullMammalRanges2, file = "~/LargeFiles/MammalRanges2.Rdata")
  
}


# Converting these to meaningful values ####

if(file.exists("~/LargeFiles/FullRangedf.Rdata")){load("~/LargeFiles/FullRangedf.Rdata"); print("Loaded!") }else{
  
  print("Converting to values!")
  
  if(file.exists("~/LargeFiles/FullValuedf1.Rdata")) load("~/LargeFiles/FullValuedf1.Rdata") else{
    
    FullValuedf <- data.frame(getValues(FullMammalRanges))
    
    print("Saving!")
    
    save(FullValuedf, file = "~/LargeFiles/FullValuedf1.Rdata")
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
      
    })
    
    save(OverList1, file = "~/LargeFiles/OverList1.Rdata")
  }
  
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
      
    })
    
    save(OverList2, file = "~/LargeFiles/OverList2.Rdata")
    
  }
  
}

print("Connecting grid dfs!")

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


Pass1Sp <- FullRangedf %>% filter(Pass == 1)
Pass2Sp <- FullRangedf %>% filter(Pass == 2)

SecondSp <- setdiff(Pass2Sp$Species, Pass1Sp$Species)

print("Making Polygons!")

if(file.exists("data/FullPolygons.Rdata")) load("data/FullPolygons.Rdata") else {
  
  FullPolygons <- lapply(names(FullMammalRanges), function(x) {
    
    print(x)
    
    r <- FullMammalRanges[[x]] > -Inf
    
    if(!all(is.na(freq(r)[,1]))) r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% mutate(Host = x) %>% return
    
  }) %>% bind_rows() %>% mutate(Pass = 1)
  
  FullPolygons2 <- lapply(SecondSp, function(x) {
    
    print(x)
    
    r <- FullMammalRanges2[[x]] > -Inf
    
    if(!all(is.na(freq(r)[,1]))) r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% mutate(Host = x) %>% return
    
  }) %>% bind_rows() %>% mutate(Pass = 2)
  
  FullPolygons <- bind_rows(FullPolygons, FullPolygons2)
  
  save(FullPolygons, file = "data/FullPolygons.Rdata")
  
}

print("Getting Range Overlap!")

if(file.exists("data/FullRangeOverlap.Rdata")) load("data/FullRangeOverlap.Rdata") else{
  
  EXT <- extent(FullMammalRanges2)
  
  FullMammalRangesb <- setExtent(FullMammalRanges, EXT, keepres = T)
  FullMammalRanges2b <- raster::subset(FullMammalRanges2, which(names(FullMammalRanges2)%in%SecondSp))
  
  MammalStack <- FullMammalRangesb
  
  for(x in names(FullMammalRanges2b)){
    
    print(x)  
    MammalStack <- raster::addLayer(MammalStack, FullMammalRanges2b[[x]])
    
  }
  
  FullRangeAdj <- PairsWisely(MammalStack)
  
  save(FullRangeAdj, file = "data/FullRangeOverlap.Rdata")
  
}

save(MammalStack, file = "data/MammalStack.Rdata")

remove(FullMammalRanges, FullMammalRanges2, MammalStack)

detach(package:raster)
