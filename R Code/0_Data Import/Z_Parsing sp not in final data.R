# Parsing spp which aren't in all mammal dataset

# Hosts in the viral data that aren't in the full mammmals #####

AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

AbsentHosts %in% Range0 # NONE of them don't have spatial data of 0 grids

AbsentHosts %in% levels(FullRangedf$Host) # ALL of them are in the range df

# Must be phylo ####

AbsentHosts %in% rownames(FullSTMatrix) # NONE with phylo data 

# Check synonyms

Synonyms <- Hosts[Hosts$Sp%in%AbsentHosts,c("Sp","Synonyms")] %>% slice(order(Sp))
Synonyms2 <- c(Synonyms$Synonyms)

names(Synonyms2) <- Synonyms$Sp

# Tadarida condylura = Mops condylurus
"Mops_condylurus" %in% rownames(FullSTMatrix)
"Micaelamys_namaquensis" %in% rownames(FullSTMatrix)
"Bos_frontalis" %in% rownames(FullSTMatrix)
"Bos_grunniens"%in% rownames(FullSTMatrix)
AbsentHosts[5] # Absent
"Capra_hircus" %in% rownames(FullSTMatrix)      # Domestic goat !!!
"Hexaprotodon_liberiensis" %in% rownames(FullSTMatrix)
"Equus_burchellii" %in% rownames(FullSTMatrix)
"Oryzomys_alfaroi"  %in% rownames(FullSTMatrix)
"Oryzomys_laticeps" %in% rownames(FullSTMatrix)
"Oryzomys_megacephalus" %in% rownames(FullSTMatrix)
"Callithrix_argentata" %in% rownames(FullSTMatrix)
"Miniopterus_schreibersii"%in% rownames(FullSTMatrix) # New species discovered; keep this one!
"Myotis_ricketti" %in% rownames(FullSTMatrix)
"Oryzomys_albigularis" %in% rownames(FullSTMatrix)
"Ovis_aries" %in% rownames(FullSTMatrix) # Sheep!
"Piliocolobus_badius"%in%rownames(FullSTMatrix)
"Piliocolobus_rufomitratus" %in%rownames(FullSTMatrix)
"Lycalopex_gymnocercus"  %in%rownames(FullSTMatrix)
"Rhinolophus_hildebrandtii"%in%rownames(FullSTMatrix)
"Oryzomys_angouya"%in%rownames(FullSTMatrix)
"Chaerephon_plicatus"%in%rownames(FullSTMatrix)
"Chaerephon_pumilus"%in%rownames(FullSTMatrix)
"Taurotragus_oryx"%in%rownames(FullSTMatrix)

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

cbind(AbsentHosts, NameReplace)

# Parsing hosts with no range size ####

H1Range0
