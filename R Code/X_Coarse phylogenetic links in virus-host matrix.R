# Coarse phylogenetic links in virus-host matrix

M %>% lapply()

vFamilyList <- 
  lapply(unique(Viruses$vFamily), function(a){ 
    lapply(unique(Hosts$hFamily), function(b) M[Viruses$vFamily == a, Hosts$hFamily==b])})


