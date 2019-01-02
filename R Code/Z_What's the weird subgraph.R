# Pruning EHA Tree ####

library(magrittr); library(igraph)

WeirdSp <- Hosts[Hosts$Degree==138, "Sp"]

Hosts[Hosts$Sp%in%WeirdSp,"hWildDomFAO"] # only one domestic
Hosts[Hosts$Sp%in%WeirdSp,"hOrder"] # Mainly bats

Hosts %>%
  subset(Sp%in%WeirdSp, select = "hOrder") %>%
  table()

Hosts %>%
  subset(!Sp%in%WeirdSp, select = "hOrder") %>%
  table()/(dim(Hosts)[1]-length(WeirdSp))*100

Hosts %>%
  subset(Sp%in%WeirdSp, select = "hOrder") %>%
  table()/(length(WeirdSp))*100

SubM <- M[,colnames(M)%in%WeirdSp]

SubM <- SubM[!rowSums(SubM)==0,]

Weirdgraph <- graph.incidence(SubM, weighted = T) 
plot(Weirdgraph) # Reckon it's rabies virus
plot(Weirdgraph, layout = layout.bipartite) # Reckon it's rabies virus

RabiesHosts <- colnames(M[,M["Rabies_virus",]>0])



