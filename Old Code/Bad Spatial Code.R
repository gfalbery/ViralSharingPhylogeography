# Spatial snippets that didn't work ####


# None of this stuff worked ####

st_crs(mammal_shapes)
class(mammal_shapes)

mammal_shapes2 <- mammal_shapes[1:100,]

plot(mammal_shapes2[,1])

i = st_intersection(mammal_shapes2[1:2,"geometry"])
i = st_difference(mammal_shapes2[1:2,"geometry"])

terr = shapefile(P("maps/Mammals_Terrestrial/TERRESTRIAL_MAMMALS.shp"), verbose = T)
terr@data$binomial = str_replace(terr@data$binomial, " ", "_")
terr = subset(terr, presence == 1)

hp3 = read_csv(P('data/hosts.csv')) %>%
  filter(hWildDomFAO == 'wild', 
         hMarOTerr == 'Terrestrial', 
         hHostNameFinal %in% spatial_host@data$binomial) %>%
  dplyr::select(hHostNameFinal) %>%
  mutate(hp3 = 1)

spatial_host <- terr

# Join spatial polygons and hp3 data frames
spatial_host@data = left_join(spatial_host@data, hp3, by = c("binomial" = 'hHostNameFinal')) %>%
  subset(hp3 == 1)

#

mammals2 = fasterize(mammal_shapes, mammal_raster, fun = "first")
mammals3 = fasterize(mammal_shapes, mammal_raster, fun = "last")
mammals4 = fasterize(mammal_shapes, mammal_raster, fun = "any")

mammals3 = fasterize(mammal_shapes, mammal_raster, field = "binomial", fun = "last")

df <- ggplot2::fortify(mammal_shapes)

#

Vertexlist <- lapply(1:2, function(x){ data.frame(fortify(mammal_shapes[x,1]), 
                                                  Sp = mammal_shapes$binomial[x],
                                                  object = mammal_shapes$objectid[x])
})

Vertexlist <- rbind(Vertexlist)

plot(df[df$binomial=="Canis lupus",1], max.plot = 26)

#

pols <- spatial_host

r <- raster(pols, res = 1)
r <- fasterize(pols, r)

# 

spatial_host@polygons[[1]]@Polygons[[2]]@coords

lin <- as(spatial_host, "SpatialLinesDataFrame")  

mpts <- raster::geom(spatial_host)

dpts <- ggplot2::fortify(spatial_host)

long_and_lat <- spatial_host %>% fortify() %>% select(long,lat)

dpts %>% 
  filter(id %in%0:100) %>% 
  ggplot(aes(long, lat, fill = id, group = group)) + 
  geom_polygon() + 
  coord_fixed() + 
  theme(legend.position = "none")