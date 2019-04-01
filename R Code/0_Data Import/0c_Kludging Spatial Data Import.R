# 0_Kludging spatial data import ####

# This involves a workaround with the SpRanger package which basically converts the rasters to a df of 1's and 0's.
# range overlap is then calculated based on joint inhabiting of these grid squares.

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

if(file.exists("data/MammalRanges.Rdata")) load(file = "data/MammalRanges.Rdata") else{
  
  mammal_shapes <- st_read("~/ShapeFiles")
  
  mammal_shapes <- st_transform(mammal_shapes, 
                                "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
  
  mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
  mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
  mammal_shapes_red <- mammal_shapes[mammal_shapes$binomial%in%unique(Hosts$Sp),]
  mammal_raster <- raster(mammal_shapes_red, res = 25000) # NB units differ from Mercator!
  
  MammalRanges <- fasterize(mammal_shapes_red, mammal_raster, by = "binomial")
  save(MammalRanges, file = "data/MammalRanges.Rdata")
  
  remove("mammal_raster", "mammal_shapes", "mammal_shapes_red")
  
}

LeftOut <- Hosts[!Hosts$Sp%in%names(MammalRanges),"Sp"]

if(file.exists("Output Files/WorldMap.Rdata")) load("Output Files/WorldMap.Rdata") else{
  
  data("wrld_simpl")
  WorldMap <- spTransform(wrld_simpl, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
  
  
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
  
}

if(!file.exists("data/RangeOverlap.Rdata")){
  
  RangeAdj <- PairsWisely(MammalRanges)
  
  save(RangeAdj, file = "data/RangeOverlap.Rdata")
  
  remove(MammalRanges)
  
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
  
  save(HostPolygons, file = "data/HostPolygons.Rdata")
  
}

# Deriving Host Centroids and Merging ####

HostCentroids <- data.frame(LongMean = with(HostPolygons, tapply(long, Host, mean)),
                            LatMean = with(HostPolygons, tapply(lat, Host, mean)),
                            Host = unique(HostPolygons$Host))

Hosts <- merge(Hosts, HostCentroids, all.x = T, by.x = "Sp", by.y = "Host")

# Viral spatial data adapted from spatial data for their hosts ####

# Making Viral Associations and Polygons ####

VirusAssocs <- apply(M, 1, function(a) names(a[a>0]))

library(centiserve); library(tidyverse); library(GGally); library(igraph)

human <- which(rownames(tFullSTMatrix) == "Homo_sapiens")

nhCytBMatrix <- tFullSTMatrix[-human,-human] # Want this to be non-human host range 

VirPhyloHostRangeMax <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(nhCytBMatrix)])>0){
    
    b <- a[a%in%rownames(nhCytBMatrix)]
    
    if(length(b)>1){
      max(nhCytBMatrix[b, b][upper.tri(nhCytBMatrix[b, b])])
    } else 0
  } else NA
  
})

VirPhyloHostRangeMean <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(nhCytBMatrix)])>0){
    
    b <- a[a%in%rownames(nhCytBMatrix)]
    
    if(length(b)>1){
      mean(nhCytBMatrix[b, b][upper.tri(nhCytBMatrix[b, b])])
    } else 0
  } else NA
  
})

VirPhyloHostRangeMedian <- sapply(VirusAssocs, function(a){
  
  if(length(a[a%in%rownames(nhCytBMatrix)])>0){
    
    b <- a[a%in%rownames(nhCytBMatrix)]
    
    if(length(b)>1){
      median(nhCytBMatrix[b, b][upper.tri(nhCytBMatrix[b, b])])
    } else 0
  } else NA
  
})

VirusHostRanges = data.frame(
  Virus = names(VirusAssocs),
  HostRangeMean = VirPhyloHostRangeMean,
  HostRangeMax = VirPhyloHostRangeMax,
  HostRangeMedian = VirPhyloHostRangeMedian
)


load("~/Albersnet/data/FullPolygons.Rdata")

CodeRoot <- "R Code/0_Data Import"

if(file.exists("~/Albersnet/data/FullRangeOverlap.Rdata")) load("~/Albersnet/data/FullRangeOverlap.Rdata") else{
  source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
}

detach(package:raster)
