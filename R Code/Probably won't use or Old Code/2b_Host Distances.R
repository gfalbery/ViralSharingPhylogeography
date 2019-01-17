
# Host Distances ####

library(igraph); library(ggregplot)

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

# Making bipartite projections ####

mFull <- table(AssocsBase[,1:2])
MFull <- as.matrix(mFull)

bipGraphFull <- graph.incidence(MFull, weighted = T)

VirusGraphFull <- bipartite.projection(bipGraphFull)$proj1
HostGraphFull <- bipartite.projection(bipGraphFull)$proj2

hHumanDist <- distances(bipGraphFull, v = V(bipGraphFull)[(dim(Viruses)[1]+1):length(V(bipGraphFull))],
                        to = V(bipGraphFull)[which(names(V(bipGraphFull)) == "Homo_sapiens")])

hHumanDist[,1] <- ifelse(hHumanDist[,1] == Inf, Inf, (hHumanDist[,1])/2)

hHumanDist2 <- data.frame(Sp = rownames(hHumanDist),
                          HumanDist = hHumanDist[,1])

hDomDist <- distances(bipGraphFull, v = V(bipGraphFull)[(dim(Viruses)[1]+1):length(V(bipGraphFull))],
                      to = V(bipGraphFull)[which(names(V(bipGraphFull)) %in% Domestics)])

hDomDist <- (hDomDist)/2

hDomDist2 <- data.frame(Sp = rownames(hDomDist),
                        MeanDomDist = apply(hDomDist, 1, function(a) mean(a[!a=="Inf"], na.rm = T)),
                        MinDomDist = apply(hDomDist, 1, function(a) ifelse(all(a==Inf),Inf,min(a[!a==Inf], na.rm = T))),
                        MaxDomDist = apply(hDomDist, 1, function(a) max(a, na.rm = T)))

Hosts <- merge(Hosts, hHumanDist2[,c("Sp","HumanDist")], by = "Sp", all.x = T)
Hosts <- merge(Hosts, hDomDist2[,c("Sp", "MeanDomDist", "MinDomDist", "MaxDomDist")], by = "Sp", all.x = T)
