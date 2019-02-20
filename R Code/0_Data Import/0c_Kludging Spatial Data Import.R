# 0_Kludging spatial data import ####

# This was all done in 2 days and is probably circuitous
# Because I mashed everything into pixels and data frames rather than working
# with spatial dataframes :howdy: :grimace:

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools)

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

CountryCentroids <- WorldMap %>% group_by(Country) %>% 
  summarise(long = mean(long), lat = mean(lat)) %>% data.frame

# THIS DATA FRAME TAKES A LOT OF MEMORY - convert to sparse matrix 
# Or learn more raster methods before pub ####

Valuedf <- data.frame(getValues(MammalRanges))
Valuedf2 <- reshape2::melt(Valuedf)
Valuedf2$x <- rep(1:MammalRanges[[1]]@ncols, MammalRanges[[1]]@nrows)
Valuedf2$y <- rep(MammalRanges[[1]]@nrows:1, each = MammalRanges[[1]]@ncols)

Valuedf3 <- data.frame(getValues(MammalRanges2))
Valuedf4 <- reshape2::melt(Valuedf3)
Valuedf4$x <- rep(1:MammalRanges2[[1]]@ncols, MammalRanges2[[1]]@nrows)
Valuedf4$y <- rep(MammalRanges2[[1]]@nrows:1, each = MammalRanges2[[1]]@ncols)

Rangedf <- rbind(Valuedf2[!is.na(Valuedf2$value),],Valuedf4[!is.na(Valuedf4$value),]) # This is where a load of them were lost ####
Rangedf <- Rangedf %>% 
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

Hosts$GeogRange <- diag(RangeOverlap)[as.character(Hosts$Sp)]

RangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)), nrow(RangeOverlap))
RangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)), nrow(RangeOverlap))

RangeAdj1 <- RangeOverlap/(RangeA + RangeB - RangeOverlap) # Weighted evenly
RangeAdj2 <- RangeOverlap/(RangeA) # Asymmetrical

Hosts[Hosts$Sp%in%rownames(RangeAdj1),"S.Greg1"] <- rowSums(RangeAdj1[as.character(Hosts[Hosts$Sp%in%rownames(RangeAdj1),"Sp"]),])
Hosts[Hosts$Sp%in%rownames(RangeAdj2),"S.Greg2"] <- rowSums(RangeAdj2[as.character(Hosts[Hosts$Sp%in%rownames(RangeAdj2),"Sp"]),])

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

if(file.exists("data/VirusPolygons.Rdata")) load("data/VirusPolygons.Rdata") else{
  
  VirusRanges <- lapply(1:length(VirusAssocs), function(b) data.frame(HostPolygons[HostPolygons$Host%in%VirusAssocs[[b]],]) %>%
                          mutate(Virus = names(VirusAssocs)[b])) %>% bind_rows()
  
  VirusPolygons <- lapply(1:length(VirusAssocs), function(x) {
    
    templist <- list()
    
    NumAssocs <- which(names(MammalRanges)%in%VirusAssocs[[x]])
    
    return(
      
      if(length(NumAssocs)>0){
        if(length(NumAssocs)>1){
          for(y in 1:length(NumAssocs)){
            r <- MammalRanges[[NumAssocs[y]]] > -Inf
            templist[[y]] <- r
          }
          
          m <- do.call(merge, templist)
          m %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
            mutate(Virus = names(VirusAssocs)[x])
          
        } else {
          r <- MammalRanges[[NumAssocs]] > -Inf
          r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
            mutate(Virus = names(VirusAssocs)[x])
        }
      } else NULL
    )
  }) 
  
  VirusPolygons <- VirusPolygons[!is.na(VirusPolygons)] %>% bind_rows()
  
  save(VirusPolygons, file = "data/VirusPolygons.Rdata")}

VirusCentroids <- data.frame(LongMean = with(VirusPolygons, tapply(long, Virus, mean)),
                             LatMean = with(VirusPolygons, tapply(lat, Virus, mean)),
                             Virus = unique(VirusPolygons$Virus))

Viruses <- merge(Viruses, VirusCentroids, all.x = T, by.x = "Sp", by.y = "Virus")

# Creating virus range overlap matrix ####

VirusRangeOverlap <- matrix(NA, nrow = nunique(Viruses$Sp), 
                            ncol = nunique(Viruses$Sp))

dimnames(VirusRangeOverlap) <- list(unique(Viruses$Sp),
                                    unique(Viruses$Sp))


if(file.exists("data/VirusRangeOverlap.Rdata")) load("data/VirusRangeOverlap.Rdata") else{
  
  GridList <- lapply(VirusAssocs, function(x){
    
    unique(Rangedf[Rangedf$Host %in% x,"GridID"]) %>% return
    
  })
  
  for(x in unique(Viruses$Sp)){
    
    VirusRangeOverlap[x,] <- sapply(GridList[unique(Viruses$Sp)], function(y) length(unlist(intersect(GridList[unique(Viruses$Sp)][[which(unique(Viruses$Sp)==x)]], y))))
    
    print(x)
    
  }
  
  diag(VirusRangeOverlap) <- sapply(GridList[unique(Viruses$Sp)], length)
  
  all(apply(VirusRangeOverlap,1,max) == diag(VirusRangeOverlap))
  
  save(GridList, file = "data/GridList.Rdata")
  save(VirusRangeOverlap, file = "data/VirusRangeOverlap.Rdata")}

VirusRangeA = matrix(rep(diag(VirusRangeOverlap), nrow(VirusRangeOverlap)), nrow(VirusRangeOverlap))
VirusRangeB = matrix(rep(diag(VirusRangeOverlap), each = nrow(VirusRangeOverlap)), nrow(VirusRangeOverlap))

VirusRangeAdj1 <- VirusRangeOverlap/(VirusRangeA + VirusRangeB - VirusRangeOverlap) # Weighted evenly
VirusRangeAdj2 <- VirusRangeOverlap/(VirusRangeA) # Asymmetrical

#load("~/Albersnet/data/FullMammalRanges.Rdata")
#load("~/Albersnet/data/FullMammalRanges2.Rdata")
load("~/Albersnet/data/FullPolygons.Rdata")

CodeRoot <- "R Code/0_Data Import"

if(file.exists("~/Albersnet/data/FullRangeOverlap.Rdata")) load("~/Albersnet/data/FullRangeOverlap.Rdata") else{
  source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
}

detach(package:raster)
