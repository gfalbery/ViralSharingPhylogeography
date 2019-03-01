

data("wrld_simpl")
WorldMap <- spTransform(wrld_simpl, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

WorldRaster <- raster(WorldMap, res = 50000) # NB units differ from Mercator!

CountryRanges <- rasterize(WorldMap, WorldRaster, by = "NAME")

CountryNames <- WorldMap@data %>% as.data.frame() %>%
  #slice(order("NAME")) %>% 
  mutate(Country = as.numeric(as.factor(NAME)))

Valuedf <- data.frame(getValues(CountryRanges))
Valuedf2 <- reshape2::melt(Valuedf)
Valuedf2$x <- rep(1:CountryRanges[[1]]@ncols, CountryRanges[[1]]@nrows)
Valuedf2$y <- rep(CountryRanges[[1]]@nrows:1, each = CountryRanges[[1]]@ncols)
Valuedf2$Country <- rep(CountryNames$Country, each = nrow(Valuedf2)/nrow(CountryNames))
Valuedf2 <- Valuedf2 %>% slice(-which(is.na(value)))

Countrydf <- Valuedf2 %>% left_join(CountryNames, by = c("value" = "Country")) %>%
  dplyr::rename(Country = NAME)

Countrydf %>% filter(Country %in%c("Italy", "Croatia")) %>%
  ggplot(aes(x, y)) + geom_tile(aes(fill = Country)) +
  coord_fixed()
