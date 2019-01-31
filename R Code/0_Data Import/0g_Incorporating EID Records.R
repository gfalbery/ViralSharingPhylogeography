
# Incorporating EID ####

EIDLocation <- read.csv("data/EID/LocationInteractions_EID2.csv", header = T)
EIDSpecies <- read.csv("data/EID/SpeciesInteractions_EID2.csv", header = T)

EveryLower <- function(a){
  A <- a %>% str_split(" ") %>% unlist
  substr(A,1,1) <- toupper(substr(A,1,1))
  return(paste(A, collapse = " "))
}

EIDLocation <- EIDLocation %>% filter(Sequences.count>0)
EIDSpecies <- EIDSpecies %>% filter(Sequences.count>0)

EIDLocation$Country <- sapply(EIDLocation$Country, EveryLower)
WorldMap2 <- WorldMap %>% group_by(Country) %>%summarise(LongMean = mean(long), LatMean = mean(lat))
WorldMap2$Country <- sapply(WorldMap2$Country, EveryLower)

setdiff(WorldMap2$Country, EIDLocation$Country)
setdiff(EIDLocation$Country, WorldMap2$Country)

EIDLocation <- merge(EIDLocation, WorldMap2, by = "Country", all.x = T) # Could get this to subregion level

Repeats <- EIDLocation$Sequences %>% str_split(";") %>% sapply(function(a) length(a[!a==""]))
Occurences <- EIDLocation$Sequences %>% str_split(";") %>% unlist

EIDIncidence <- EIDLocation[rep(1:nrow(EIDLocation), Repeats),]
EIDIncidence$Sequence <- Occurences

# So now we have a spatiotemporal record of every sequence + country etc.

# Doing the same w species associations ####

Repeats <- EIDSpecies$Sequences %>% str_split(";") %>% sapply(function(a) length(a[!a==""]))
Occurences <- EIDSpecies$Sequences %>% str_split(";") %>% unlist

EIDSpecies2 <- EIDSpecies[rep(1:nrow(EIDSpecies), Repeats),]
EIDSpecies2$Sequence <- Occurences

EIDIncidence <- merge(EIDIncidence, EIDSpecies, by = "Sequences")

table(EIDIncidence$Cargo==as.character(EIDIncidence$Species)) # uh oh
EIDIncidence[which(EIDIncidence$Cargo!=as.character(EIDIncidence$Species)),c("Cargo","Species")]

ggplot(EIDIncidence, aes(LongMean, LatMean)) + 
  geom_point(alpha = 0.3, position = position_jitter(w = 100000, h = 100000)) + 
  coord_fixed() + facet_wrap(~Cargo.classification)

EIDViruses <- EIDIncidence %>% filter(Cargo.classification=="Virus")
dim(EIDViruses)
M <- table(EIDViruses$Carrier,EIDViruses$Cargo)

table(c(M)>0)


