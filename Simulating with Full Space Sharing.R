
# Simulating if all species overlapped ####

# Simulating on the full network ####

library(MCMCglmm); library(tidyverse); library(Matrix); library(parallel)

# Trying it without random effects ####

AllSpacePredList <- list()

MX1s = model.matrix( ~ Space + Phylo + Space:Phylo, data = mutate(AllMammaldf, Space=max(Space)) ) %>%
  as.matrix %>% as("dgCMatrix")

N = nrow(AllMammaldf)

ClusterMCMC <- p$df

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

AllSpacePredList <- mclapply(1:length(RowsSampled), function(x){
  
  RowSampled <- RowsSampled[x]
  
  XFX <- p$df[RowSampled, XBetas] %>% unlist
  XPredictions <- c(XFX %*% t(MX1s))
  
  ZPredictions2a <- rnorm(n = N, mean = mean(d$DCites)*p$df[RowSampled,"beta_d_cites_s"], sd = p$df[RowsSampled[x], "sigma"])
  ZPredictions2b <- rnorm(n = N, mean = mean(d$DCites)*p$df[RowSampled,"beta_d_cites_s"], sd = p$df[RowsSampled[x], "sigma"])
  
  Predictions <- XPredictions[[1]]@x + ZPredictions2a + ZPredictions2b
  
  PZero <- rbinom(n = N, size = 1, prob = logistic(Predictions))
  
  PZero
  
}, mc.cores = 10)

#AllSpacePredList <- AllSpacePredList[-which(sapply(AllSpacePredList, function(a) is.null(a)|class(a)=="try-error")]

AllSpacePredDF <- as.data.frame(AllSpacePredList) #%>% as.matrix %>% unname %>% as("dgCMatrix")
AllMammaldf$SpacePredVirus <- apply(AllSpacePredDF,1, function(a) a %>% mean)

# Simulating the network #####

AllSpaceSims <- parallel::mclapply(1:length(AllSpacePredList), function(i){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperMammals)] <- round(AllSpacePredList[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                           union(AllMammaldf$Sp,AllMammaldf$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

AllSpaceSimGs <- parallel::mclapply(1:length(AllSpacePredList), function(i){
  
  graph.adjacency(AllSpaceSims[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

AllSpaceSimGs <- AllSpaceSimGs[-which(sapply(AllSpaceSims, is.null))]

AllSpaceDegdf <- sapply(AllSpaceSimGs, function(a) degree(a)) #%>% as.data.frame
AllSpaceEigendf <- sapply(AllSpaceSimGs, function(a) eigen_centrality(a)$vector)# %>% as.data.frame

AllSpacePredDegrees <- apply(AllSpaceDegdf, 1, mean)
AllPredEigen <- apply(AllEigendf, 1, mean)

AllSpaceComponents <- sapply(AllSpaceSimGs, function(b) components(b)$no)
AllSpaceCluster = sapply(AllSpaceSimGs, function(a) transitivity(a)) # all zero, don't bother

AllBetweenness = parallel::mclapply(AllSimGs, function(a) betweenness(a), mc.cores = 10)

AllCloseness = parallel::mclapply(AllSimGs, function(a) closeness(a), mc.cores = 10)

Hosts[,"AllSpacePredDegree"] <- AllSpacePredDegrees[as.character(Hosts$Sp)]

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")
Panth1[,"AllSpacePredDegree"] <- AllSpacePredDegrees[as.character(Panth1$Sp)]

#FullPolygons[,"AllPredDegree"] <- AllPredDegrees[as.character(FullPolygons$Host)]

#ggplot(FullPolygons, aes(long, lat, group = paste(group,Host))) + geom_polygon(aes(fill = AllPredDegree), alpha = 0.01)  + 
#  coord_fixed() + 
#  theme_void() + 
#  scale_fill_gradient(low = "white", high = AlberColours[3]) +  
#  theme(legend.position = "none") +
#  ggsave("Figures/Degree Averaging.jpeg", units = "mm", height = 100, width = 100, dpi = 300)

# Identifying within-order links ####

Panth1 <- Panth1 %>% slice(order(Panth1$Sp)) %>% 
  filter(Sp%in%V(AllSimGs[[1]])$name) %>% droplevels

#all(Panth1$Sp==V(AllSimGs[[1]])$name)

hOrderList <- list()

hOrderList <- mclapply(1:length(AllSimGs), function(i){
  lapply(levels(Panth1$hOrder), function(x) induced_subgraph(AllSimGs[[i]], 
                                                             Panth1$hOrder == x))
})


for(x in levels(Panth1$hOrder)){
  hOrderList[[x]] <- induced_subgraph(AllSimGs[[1]], 
                                      Panth1$hOrder == x)
}






