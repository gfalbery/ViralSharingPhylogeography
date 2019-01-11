# 0_Kludging spatial data import ####

# This was all done in 2 days and is probably circuitous
# Because I mashed everything into pixels and data frames rather than working
# with spatial dataframes :howdy: :grimace:

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools)

if(file.exists("data/MammalRanges.Rdata")) load(file = "data/MammalRanges.Rdata") else{
  
  mammal_shapes <- st_read("maps/Mammals_Terrestrial")
  
  mammal_shapes <- st_transform(mammal_shapes, 54009) # Mollweide projection 
  
  # Mollweide projection = +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
  # This projection retains grid size as much as possible, but at the expense of shape
  
  mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
  mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
  mammal_shapes_red <- mammal_shapes[mammal_shapes$binomial%in%unique(Hosts$Sp),]
  mammal_raster <- raster(mammal_shapes_red, res = 50000) # NB units differ from Mercator!
  
  MammalRanges <- fasterize(mammal_shapes_red, mammal_raster, by = "binomial")
  #save(MammalRanges, file = "data/MammalRanges.Rdata")
  
  AllMammals <- fasterize(mammal_shapes_red[-which(mammal_shapes_red$binomial == "Ursus_maritimus"),], mammal_raster, fun = "sum")
  
}

data("wrld_simpl")
WorldMap <- spTransform(wrld_simpl, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
#WorldMapRaster <- raster(WorldMap, res = 50000) # NB units differ from Mercator!
#World <- fasterize()

# I HATE this 

WorldPolygons <- lapply(WorldMap@polygons, 
                        function(b){ lapply(b@Polygons, 
                                            function(a) as.data.frame(a@coords))}) %>% 
  unlist(recursive = F)

for(x in 1:length(WorldPolygons)) WorldPolygons[[x]]$group <- x 

WorldMap <- bind_rows(WorldPolygons) %>% rename(long = V1, lat = V2)

ggplot(WorldMap, aes(long, lat, group = group)) + 
  geom_polygon(colour = "dark grey", fill = NA) +
  coord_fixed()

remove("mammal_raster", "mammal_shapes", "mammal_shapes_red")

# THIS DATA FRAME TAKES A LOT OF MEMORY - convert to sparse matrix 
# Or learn more raster methods before pub ####

Valuedf <- data.frame(getValues(MammalRanges))
Valuedf2 <- reshape2::melt(Valuedf)
Valuedf2$x <- rep(1:MammalRanges[[1]]@ncols, MammalRanges[[1]]@nrows)
Valuedf2$y <- rep(MammalRanges[[1]]@nrows:1, each = MammalRanges[[1]]@ncols)
Rangedf <- Valuedf2[!is.na(Valuedf2$value),]
Rangedf <- Rangedf %>% 
  rename(Host = variable, Presence = value)

Rangedf %>% filter(Host %in% levels(Hosts$Sp)) %>% ggplot(aes(x, y, fill = Host)) + geom_tile() + 
  coord_fixed() +
  lims(x = c(1, MammalRanges[[1]]@ncols), y = c(1, MammalRanges[[1]]@nrows)) +
  theme(legend.position = "none")

Rangedf$GridID <- with(Rangedf, paste(x, y))

Range0 <- levels(Rangedf$Host)[which(table(Rangedf$Host)==0)] # Hosts that have no spatial records??
Rangedf <- droplevels(Rangedf) 

# Using igraph to project it

# M2 <- with(Rangedf, table(Host, GridID))
# spacebipgraph <- graph.incidence(M2, weighted = T)
# SpaceHostGraph <- bipartite.projection(spacebipgraph)$proj2
# SpaceHostadj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))

# Making home range overlaps ####
# There MUST be a quicker way of doing this but for now I don't have it

RangeOverlap <- matrix(0, nrow = nlevels(Rangedf$Host), ncol = nlevels(Rangedf$Host))
dimnames(RangeOverlap) <- list(levels(Rangedf$Host),levels(Rangedf$Host))

for(x in levels(Rangedf$Host)){
  
  Grids <- Rangedf[Rangedf$Host==x,"GridID"]
  SubRangedf <- Rangedf[Rangedf$GridID %in% Grids,]
  
  RangeOverlap[x,] <- table(SubRangedf$Host)
  
  print(x)
  
}

diag(RangeOverlap) # Range Size for each species

RangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)), nrow(RangeOverlap))
RangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)), nrow(RangeOverlap))

RangeAdj1 <- RangeOverlap/(RangeA + RangeB - RangeOverlap) # Weighted evenly
RangeAdj2 <- RangeOverlap/(RangeA) # Asymmetrical

# Making polygons for display ####

HostPolygons <- lapply(levels(Valuedf2$variable), function(x) {
  
  if(!x%in%Range0){
    
    r <- MammalRanges[[x]] > -Inf
    
    r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
      mutate(Host = x) %>% return
  }
  
})

HostPolygons <- bind_rows(HostPolygons)

ggplot(HostPolygons, aes(long, lat, group = group)) + 
  geom_path(data = WorldMap) +
  geom_path(alpha = 0.6, aes(colour = Host)) + coord_fixed() + theme(legend.position = "none")

# Deriving centroids ####

HostCentroids <- data.frame(LongMean = with(HostPolygons, tapply(long, Host, mean)),
                            LatMean = with(HostPolygons, tapply(lat, Host, mean)),
                            Host = unique(HostPolygons$Host))

# Viral spatial data ####

VirusAssocs <- apply(M, 1, function(a) names(a[a>0]))

VirusRanges <- lapply(1:length(VirusAssocs), function(b) data.frame(HostPolygons[HostPolygons$Host%in%VirusAssocs[[b]],]) %>%
                        mutate(Virus = names(VirusAssocs)[b])) %>% bind_rows()

VirusPolygons <- lapply(1:length(VirusAssocs), function(x) {
  
  templist <- list()
  
  NumAssocs <- which(unique(mammal_shapes_red$binomial)%in%VirusAssocs[[x]])
  
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

ggplot(VirusPolygons[VirusPolygons$Virus%in%head(unique(VirusPolygons$Virus), 25),], 
       aes(long, lat, colour = Virus, group = paste(Virus, group))) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_path() + 
  facet_wrap(~Virus) + coord_fixed() +
  ggtitle("Geographic Ranges of Viral Hosts") + theme(legend.position = "none") +
  ggsave("Figures/Viral Host Ranges.jpeg", units = "mm", width = 250, height = 150, dpi = 300)

VirusCentroids <- data.frame(LongMean = with(VirusPolygons, tapply(long, Virus, mean)),
                             LatMean = with(VirusPolygons, tapply(lat, Virus, mean)),
                             Virus = unique(VirusPolygons$Virus))

VirusCentroids2 <- data.frame(LongMean = with(VirusRanges, tapply(long, Virus, mean)),
                              LatMean = with(VirusRanges, tapply(lat, Virus, mean)),
                              Virus = unique(VirusRanges$Virus))

# Merging ####

SpatialHosts <- merge(Hosts, HostCentroids, by.x = "Sp", by.y = "Host")
SpatialViruses <- merge(Viruses, VirusCentroids, by.x = "Sp", by.y = "Virus")

VirusRangeOverlap <- matrix(0, nrow = nunique(SpatialViruses$Sp), 
                            ncol = nunique(SpatialViruses$Sp))

dimnames(VirusRangeOverlap) <- list(unique(SpatialViruses$Sp),
                                    unique(SpatialViruses$Sp))

GridList <- lapply(VirusAssocs, function(x){
  
  unique(Rangedf[Rangedf$Host %in% x,"GridID"]) %>% return
  
})


for(x in unique(SpatialViruses$Sp)){
  
  VirusRangeOverlap[x,] <- sapply(GridList[unique(SpatialViruses$Sp)], function(y) length(unlist(intersect(GridList[unique(SpatialViruses$Sp)][[which(unique(SpatialViruses$Sp)==x)]], y))))
  
}

diag(VirusRangeOverlap) <- sapply(GridList[unique(SpatialViruses$Sp)], length)

all(apply(VirusRangeOverlap,1,max)==diag(VirusRangeOverlap))

VirusRangeA = matrix(rep(diag(VirusRangeOverlap), nrow(VirusRangeOverlap)), nrow(VirusRangeOverlap))
VirusRangeB = matrix(rep(diag(VirusRangeOverlap), each = nrow(VirusRangeOverlap)), nrow(VirusRangeOverlap))

VirusRangeAdj1 <- VirusRangeOverlap/(VirusRangeA + VirusRangeB - VirusRangeOverlap) # Weighted evenly
VirusRangeAdj2 <- VirusRangeOverlap/(VirusRangeA) # Asymmetrical

VirusAssocs[unique(SpatialViruses$Sp)][which(sapply(GridList[unique(SpatialViruses$Sp)], length)==0)]
# Those that have no distribution data include mainly viruses of domestic, human, marine hosts

# Plotting out ####

ggplot(HostCentroids, aes(LongMean, LatMean)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(alpha = 0.6, colour = AlberColours[3]) + 
  coord_fixed() + ggtitle("Host Geographic Centroids") +
  labs(x = "Longitude", y = "Latitude") +
  ggsave("Figures/HostCentroids.jpg", units = "mm", width = 150, height = 80)

ggplot(VirusCentroids, aes(LongMean, LatMean)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(alpha = 0.6, colour = AlberColours[5]) + 
  coord_fixed() + ggtitle("Virus Geographic Centroids") +
  labs(x = "Longitude", y = "Latitude") +
  ggsave("Figures/VirusCentroids.jpg", units = "mm", width = 150, height = 80)

# I wonder if centrality in the network is spatially autocorrelated?

ggplot(SpatialHosts, aes(LongMean, LatMean)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(aes(size = Eigenvector), alpha = 0.6, colour = AlberColours[3]) + 
  coord_fixed() + ggtitle("Host Location:Centrality") +
  labs(x = "Longitude", y = "Latitude") +
  ggsave("Figures/HostCentroids2.jpg", units = "mm", width = 150, height = 80)

ggplot(SpatialViruses, aes(LongMean, LatMean)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(aes(size = vEigenvector), alpha = 0.6, colour = AlberColours[5]) + 
  coord_fixed() + ggtitle("Virus Location:Centrality") +
  labs(x = "Longitude", y = "Latitude") +
  #ggrepel::geom_text_repel(data = SpatialViruses %>% filter(vEigenvector%in%Largest(SpatialViruses$Eigenvector)),aes(label = Sp)) +
  ggsave("Figures/VirusCentroids2.jpg", units = "mm", width = 150, height = 80)

# Combining Centroids and ranges ####

ggplot(SpatialHosts, aes(LongMean, LatMean)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_path(data = HostPolygons, aes(long, lat, colour = Host, group = paste(Host, group)), alpha = 0.6) + 
  geom_point(alpha = 0.6, colour = "black") + 
  coord_fixed() + ggtitle("Host Ranges") +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  ggsave("Figures/HostCentroids3.jpg", units = "mm", width = 150, height = 80)

ggplot(SpatialViruses, aes(LongMean, LatMean)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_path(data = VirusPolygons, 
            aes(long, lat, colour = Virus, group = paste(Virus, group)), alpha = 0.6) + 
  geom_point(alpha = 0.6, colour = "black") + 
  coord_fixed() + ggtitle("Virus Ranges") +
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  ggsave("Figures/VirusCentroids3.jpg", units = "mm", width = 150, height = 80)
