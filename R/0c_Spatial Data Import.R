# Creating exhaustive mammal spatial data ####

# Rscript "R Code/Spatial.R"

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Importing/making ranges ####

if(file.exists("data/FullRangeOverlap.Rdata")){ load("data/FullRangeOverlap.Rdata"); load("data/FullPolygons.Rdata") }else{
  
  print("Mammal Ranges 1!")
  
  if(file.exists("~/LargeFiles/MammalRanges1.Rdata")) load("~/LargeFiles/MammalRanges1.Rdata") else{
    
    mammal_shapes <- st_read("~/ShapeFiles")
    
    Transform = T
    
    if(Transform){
      
      mammal_shapes <- st_transform(mammal_shapes, 
                                    "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
    } 
    
    mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
    mammal_shapes <- mammal_shapes[order(mammal_shapes$binomial),]
    mammal_raster_full <- raster(mammal_shapes, res = 0.25)
    
    print("Fasterising!")
    
    FullMammalRanges <- fasterize(mammal_shapes, mammal_raster_full, by = "binomial")
    
    if(!Transform) FullMammalRanges <- FullMammalRanges*AreaRaster
    
    names(FullMammalRanges) <- unique(mammal_shapes$binomial)
    
    save(FullMammalRanges, file = "~/LargeFiles/MammalRanges1.Rdata")
    
  }
  
  # Trying earlier dataset ####
  
  if(file.exists("~/LargeFiles/MammalRanges2.Rdata")) load("~/LargeFiles/MammalRanges2.Rdata") else{
    
    print("Mammal Ranges 2!")
    
    mammal_shapes2 <- st_read("~/ShapeFiles2")
    
    if(Transform){
      mammal_shapes2 <- st_transform(mammal_shapes2, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection
    }
    
    mammal_shapes2$binomial = str_replace(mammal_shapes2$BINOMIAL, " ", "_")
    mammal_shapes2 <- mammal_shapes2[order(mammal_shapes2$binomial),]
    
    mammal_raster_full2 <- raster(mammal_shapes2, res = 0.25) # NB units differ from Mercator!
    
    print("Fasterising!")
    
    FullMammalRanges2 <- fasterize(mammal_shapes2, mammal_raster_full2, by = "binomial")
    
    if(!Transform){
      
      FullMammalRanges2 <- FullMammalRanges2*AreaRaster
      
    }
    
    names(FullMammalRanges2) <- unique(mammal_shapes2$binomial)
    
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
  
  Pass1Sp <- names(FullMammalRanges)
  Pass2Sp <- names(FullMammalRanges2)
  
  SecondSp <- setdiff(Pass2Sp, Pass1Sp)
  
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
    
    MammalStack <- lapply(names(FullMammalRanges), function(a){
      print(a)
      FullMammalRanges[[a]] %>% resample(blank, method = 'ngb')
    })
    
    MammalStack2 <- lapply(SecondSp, function(a){
      print(a)
      FullMammalRanges2[[a]] %>% resample(blank, method = 'ngb')
    })
    
    names(MammalStack2) <- SecondSp
    
    MammalStackFull <- append(MammalStack, MammalStack2)
    names(MammalStackFull) <- c(Pass1Sp, SecondSp)
    
    MammalStackFull <- MammalStackFull[sort(names(MammalStackFull))]
    
    FullRangeAdj <- PairsWisely(MammalStackFull, Area = T)
    
    save(FullRangeAdj, file = "data/FullRangeOverlap.Rdata")
    
  }
  
  save(MammalStackFull, file = "data/MammalStack.Rdata")
  
  detach(package:raster)
  remove(FullMammalRanges, FullMammalRanges2, MammalStack)
  
}

# 0_Kludging spatial data import ####

# This involves a workaround with the SpRanger package which basically converts the rasters to a df of 1's and 0's.
# range overlap is then calculated based on joint inhabiting of these grid squares.

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

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

RangeHosts <- intersect(Hosts$Sp, rownames(FullRangeAdj))

RangeAdj <- FullRangeAdj[RangeHosts,RangeHosts]

# Viral spatial data adapted from spatial data for their hosts ####

# Making Viral Associations and Polygons ####

VirusAssocs <- apply(M, 1, function(a) names(a[a>0]))

detach(package:raster)

library(centiserve); library(tidyverse); library(GGally); library(igraph)

# Removing non-eutherian mammals ####

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

FinalHostNames <- reduce(list(
  rownames(FullRangeAdj), 
  colnames(FullSTMatrix),
  rownames(HostAdj)), intersect)

FHN <- FinalHostNames %>% setdiff(NonEutherianSp); length(FHN)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj))
AllMammals <- AllMammals[order(AllMammals)]
AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

# Replacing absent names in the full ST matrix ####

NameReplace <- c(
  "Micaelamys_namaquensis",
  "Akodon_paranaensis",
  "Bos_frontalis",
  "Bos_grunniens",
  "Bubalus_arnee", # Absent
  "Capra_hircus",
  "Hexaprotodon_liberiensis",
  "Equus_burchellii",
  "Oryzomys_alfaroi" ,
  "Oryzomys_laticeps",
  "Oryzomys_megacephalus",
  "Callithrix_argentata",
  "Miniopterus_schreibersii",
  "Myotis_ricketti",
  "Oryzomys_albigularis",
  "Ovis_aries",
  "Piliocolobus_badius",
  "Piliocolobus_rufomitratus" ,
  "Lycalopex_gymnocercus" ,
  "Rhinolophus_hildebrandtii",
  "Oryzomys_angouya",
  "Mops_condylurus",
  "Chaerephon_plicatus",
  "Chaerephon_pumilus",
  "Taurotragus_oryx")

names(NameReplace) <- AbsentHosts

rownames(FullSTMatrix) <- colnames(FullSTMatrix) <- sapply(rownames(FullSTMatrix), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

tFullSTMatrix <- 1 - (FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp] - 
                        min(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp]))/
  max(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp])

tSTMatrix <- tFullSTMatrix

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


