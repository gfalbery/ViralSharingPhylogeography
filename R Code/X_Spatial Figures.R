# X_Spatial Figures ####

# Plotting out ####

ggplot(VirusPolygons[VirusPolygons$Virus%in%head(unique(VirusPolygons$Virus), 25),], 
       aes(long, lat, colour = Virus, group = paste(Virus, group))) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_path() + 
  facet_wrap(~Virus) + coord_fixed() +
  ggtitle("Geographic Ranges of Viral Hosts") + theme(legend.position = "none") +
  ggsave("Figures/Viral Host Ranges.jpeg", units = "mm", width = 250, height = 150, dpi = 300)

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



ggplot(SpatialViruses, aes(LatMean.Total, LatMean.Centroid)) + geom_point() + coord_fixed()
ggplot(SpatialViruses, aes(LongMean.Total, LongMean.Centroid)) + geom_point() + coord_fixed()

ggplot(SpatialViruses, aes(LongMean.Centroid, LatMean.Centroid)) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(alpha = 0.6, colour = AlberColours[5]) + 
  coord_fixed() + ggtitle("Virus Spatial Centroids of Centroid Host") +
  labs(x = "Longitude", y = "Latitude")


