
# References as sampling events ####

Refs <- read.csv("D:/Reference Sources.csv", header = T)

Refs <- Refs %>% rename(Name = "ï..Name", Year = "publication_year") %>%
  select(-c(""))

setdiff(Refs$sampling_location, CountryCentroids$Country)
setdiff(CountryCentroids$Country, Refs$sampling_location)

Refs <- merge(Refs, CountryCentroids, 
              by.x = "sampling_location", by.y  = "Country", all.x = T) 

table(is.na(Refs$long))

ggplot(Refs, aes(long, lat, colour = Year)) + 
  geom_point(position = "jitter") +
  coord_fixed() +
  facet_wrap(~Year)

