
# Host Distances ####

library(igraph); library(ggregplot)

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

hHumanDist <- distances(bipgraph, v = V(bipgraph)[(dim(Viruses)[1]+1):length(V(bipgraph))],
                        to = V(bipgraph)[which(names(V(bipgraph)) == "Homo_sapiens")], weights=NA)

hHumanDist[,1] <- ifelse(hHumanDist[,1] == Inf, Inf, (hHumanDist[,1])/2)

hHumanDist2 <- data.frame(Sp = rownames(hHumanDist),
                          HumanDist = hHumanDist[,1],
                          hDom = Hosts$hDom,
                          Family = Hosts$hFamily)

hDomDist <- distances(bipgraph, v = V(bipgraph)[(dim(Viruses)[1]+1):length(V(bipgraph))],
                      to = V(bipgraph)[which(names(V(bipgraph)) %in% Domestics)], weights=NA)

hDomDist <- (hDomDist)/2

hDomDist2 <- data.frame(Sp = rownames(hDomDist),
                        MeanDomDist = apply(hDomDist, 1, function(a) mean(a[!a=="Inf"], na.rm = T)),
                        MinDomDist = apply(hDomDist, 1, function(a) ifelse(all(a==Inf),Inf,min(a[!a==Inf], na.rm = T))),
                        MaxDomDist = apply(hDomDist, 1, function(a) max(a, na.rm = T)),
                        Family = Hosts$hFamily)

Hosts <- merge(Hosts, hHumanDist2[,c("Sp","HumanDist")], by = "Sp", all.x = T)
Hosts <- merge(Hosts, hDomDist2[,c("Sp", "MeanDomDist", "MinDomDist", "MaxDomDist")], by = "Sp", all.x = T)

