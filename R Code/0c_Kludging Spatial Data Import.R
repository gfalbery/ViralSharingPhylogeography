# 0_Kludging spatial data import ####

# This was all done in 2 days and is probably circuitous
# Because I mashed everything into pixels and data frames rather than working
# With spatial dataframes :howdy:

library(sf); library(fasterize); library(Matrix);
library(ggregplot); library(raster); library(tidyverse); library(igraph)

if(file.exists("data/MammalRanges.Rdata")) load(file = "data/MammalRanges.Rdata") else{
  
  mammal_shapes <- st_read("maps/Mammals_Terrestrial")
  
  mammal_shapes <- st_transform(mammal_shapes, 54009) # Mollweide projection
  # This projection retains grid size as much as possible, but at the expense of shape
  
  mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
  mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
  mammal_shapes_red <- mammal_shapes[mammal_shapes$binomial%in%unique(Hosts$Sp),]
  mammal_raster <- raster(mammal_shapes_red, res = 50000) # NB units differ from Mercator!
  
  MammalRanges = fasterize(mammal_shapes_red, mammal_raster, by = "binomial")
  #save(MammalRanges, file = "data/MammalRanges.Rdata")
  
  AllMammals = fasterize(mammal_shapes_red[-which(mammal_shapes_red$binomial == "Ursus_maritimus"),], mammal_raster, fun = "sum")
  
}

WorldOutline <- AllMammals %>% rasterToPolygons(dissolve=TRUE) %>% fortify

ggplot(WorldOutline, aes(long, lat, group = group, colour = rownames(WorldOutline))) + 
  geom_path() + theme_void() + theme(legend.position = "none") + 
  coord_fixed() +
  ggsave("Figures/CoolMap.tiff", units = "mm", width = 250, height = 120, dpi = 300)

ggplot(WorldOutline, aes(long, lat, group = group, colour = rownames(WorldOutline))) + 
  geom_path(size = 0.25) + theme_void() + theme(legend.position = "none") + 
  coord_fixed() +
  ggsave("Figures/CoolMap.jpeg", units = "mm", width = 250, height = 120, dpi = 300)

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
RangeAdj2 <- RangeOverlap/(RangeB) # Asymmetrical

# Making polygons for display ####

HostPolygons <- lapply(levels(Valuedf2$variable), function(x) {
  
  if(!x%in%Range0){
    
    r <- MammalRanges[[x]] > -Inf
    
    r %>% rasterToPolygons(dissolve=TRUE) %>% fortify %>% 
      mutate(Host = x) %>% return
  }
  
})

HostPolygons <- bind_rows(HostPolygons)

ggplot(HostPolygons, aes(long, lat, colour = Host, group = paste(Host, group))) + 
  geom_path(alpha = 0.6) + coord_fixed() + theme(legend.position = "none")

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
  
  VirusRangeOverlap[x,] <- sapply(GridList[unique(SpatialViruses$Sp)], function(y) length(unlist(intersect(GridList[unique(SpatialViruses$Sp)][[x]], y))))
  
}

diag(VirusRangeOverlap) <- sapply(GridList[unique(SpatialViruses$Sp)], length)

# Plotting out ####

ggplot(HostCentroids, aes(LongMean, LatMean)) + 
  geom_point(alpha = 0.6, colour = AlberColours[3]) + 
  coord_fixed() + ggtitle("Host Geographic Centroids") +
  ggsave("Figures/HostCentroids.jpg", units = "mm", width = 150, height = 80)

ggplot(VirusCentroids, aes(LongMean, LatMean)) + 
  geom_point(alpha = 0.6, colour = AlberColours[5]) + 
  coord_fixed() + ggtitle("Virus Geographic Centroids") +
  ggsave("Figures/VirusCentroids.jpg", units = "mm", width = 150, height = 80)

# I wonder if centrality in the network is spatially autocorrelated?

ggplot(SpatialHosts, aes(LongMean, LatMean)) + 
  geom_point(aes(size = Eigenvector), alpha = 0.6, colour = AlberColours[3]) + 
  coord_fixed() + ggtitle("Host Location:Centrality") +
  ggsave("Figures/HostCentroids2.jpg", units = "mm", width = 150, height = 80)

ggplot(SpatialViruses, aes(LongMean, LatMean)) + 
  geom_point(aes(size = vEigenvector), alpha = 0.6, colour = AlberColours[5]) + 
  coord_fixed() + ggtitle("Virus Location:Centrality") +
  ggsave("Figures/VirusCentroids2.jpg", units = "mm", width = 150, height = 80)

# Combining Centroids and ranges ####

ggplot(SpatialHosts, aes(LongMean, LatMean)) + 
  geom_path(data = HostPolygons, aes(long, lat, colour = Host, group = paste(Host, group)), alpha = 0.6) + 
  geom_point(alpha = 0.6, colour = "black") + 
  coord_fixed() + ggtitle("Host Ranges") +
  theme(legend.position = "none") +
  ggsave("Figures/HostCentroids3.jpg", units = "mm", width = 150, height = 80)

ggplot(SpatialViruses, aes(LongMean, LatMean)) + 
  geom_path(data = VirusPolygons, 
            aes(long, lat, colour = Virus, group = paste(Virus, group)), alpha = 0.6) + 
  geom_point(alpha = 0.6, colour = "black") + 
  coord_fixed() + ggtitle("Virus Ranges") +
  theme(legend.position = "none") +
  ggsave("Figures/VirusCentroids3.jpg", units = "mm", width = 150, height = 80)
