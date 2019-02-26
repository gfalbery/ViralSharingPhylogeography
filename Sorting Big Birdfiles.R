
# Sorting Big Birdfiles

library(fs)

setwd("~")

Paths <- dir_info() %>% filter(modification_time<"2018-12-01") %>% dplyr::select(path) %>% as.list()

Destination <- dir_info() %>% filter(path == "AllBirds") %>% dplyr::select(path)

file.copy(Paths$path, Destination$path)

file.remove(Paths$path)

dir_info("AllBirds") %>% dplyr::select(path) %>% str_split("/")


# Getting all bird shapefiles in ####
getwd()

Bird_shapes <- st_read("~/AllBirds")

#Bird_shapes <- st_transform(Bird_shapes, 54009) # Mollweide projection 
Bird_shapes <- st_transform(Bird_shapes, 
                            "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 

# Mollweide projection = +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs
# This projection retains grid size as much as possible, but at the expense of shape

Bird_shapes$binomial = str_replace(Bird_shapes$binomial, " ", "_")
Bird_shapes <- Bird_shapes[order(Bird_shapes$binomial),]
Bird_shapes_red <- Bird_shapes[Bird_shapes$binomial%in%unique(Hosts$Sp),]
Bird_raster <- raster(Bird_shapes_red, res = 50000) # NB units differ from Mercator!

BirdRanges <- fasterize(Bird_shapes_red, Bird_raster, by = "binomial")
#save(BirdRanges, file = "data/BirdRanges.Rdata")