
# Simulating on the full network ####

library(MCMCglmm); library(tidyverse); library(Matrix); library(parallel)

if(file.exists("~/Albersnet/Bin Model Output.Rdata")) load("~/Albersnet/Bin Model Output.Rdata") else{
  
  p <- process_stanfit(BinModel, n.pars.to.trim = 3) # Takes a WHILE
  
}

invlogit <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj1[AllMammals,AllMammals]),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
)

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals,AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

# Trying it without random effects ####

MX1 = model.matrix( ~ Space + Phylo + Space:Phylo, data = AllMammaldf) %>%
  as.matrix %>% as("dgCMatrix")

N = nrow(AllMammaldf)

ClusterMCMC <- p$df

RowsSampled <- sample(1:nrow(ClusterMCMC), 100, replace = F)

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

AllPredList <- mclapply(1:length(RowsSampled), function(x){
  
  RowSampled <- RowsSampled[x]
  
  XFX <- p$df[RowSampled, XBetas] %>% unlist
  XPredictions <- c(XFX %*% t(MX1))
  
  ZPredictions2a <- rnorm(n = N, mean = mean(d$DCites)*p$df[RowSampled,"beta_d_cites_s"], sd = p$df[RowsSampled[x], "sigma"])
  ZPredictions2b <- rnorm(n = N, mean = mean(d$DCites)*p$df[RowSampled,"beta_d_cites_s"], sd = p$df[RowsSampled[x], "sigma"])
  
  Predictions <- XPredictions[[1]]@x + ZPredictions2a + ZPredictions2b
  
  PZero <- rbinom(n = N, size = 1, prob = logistic(Predictions))
  
  PZero
  
}, mc.cores = 10)

AllPredDF <- as.data.frame(AllPredList)
AllMammaldf$PredVirus <- apply(AllPredDF,1, function(a) a %>% mean)

AllMammaldf$PredVirusQ <- cut(AllMammaldf$PredVirus,
                              breaks = c(-1:10/10),
                              labels = c(0:10/10))

# Simulating the network #####

AllSims <- parallel::mclapply(1:length(AllPredList), function(i){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                   ncol = length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperMammals)] <- round(AllPredList[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                           union(AllMammaldf$Sp,AllMammaldf$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

AllSimGs <- parallel::mclapply(1:length(AllPredList), function(i){
  
  graph.adjacency(AllSims[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

if(length(which(sapply(AllSims, is.null)))>0) AllSimGs <- AllSimGs[-which(sapply(AllSims, is.null))]

AllPrev = sapply(AllPredList, Prev)

AllDegdf <- sapply(AllSimGs, function(a) degree(a))
AllPredDegrees <- apply(AllDegdf, 1, mean)

#AllEigendf <- sapply(AllSimGs, function(a) eigen_centrality(a)$vector)
#AllPredEigen <- apply(AllEigendf, 1, mean)

AllComponents <- sapply(AllSimGs, function(b) components(b)$no)

#AllCluster = sapply(AllSimGs, function(a) transitivity(a))
#Betweenness = parallel::mclapply(AllSimGs, function(a) betweenness(a), mc.cores = 10)
#AllCloseness = parallel::mclapply(AllSimGs, function(a) closeness(a), mc.cores = 10)

Hosts[,"AllPredDegree"] <- AllPredDegrees[as.character(Hosts$Sp)]

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")
Panth1[,"AllPredDegree"] <- AllPredDegrees[as.character(Panth1$Sp)]

FullPolygons[,"AllPredDegree"] <- AllPredDegrees[as.character(FullPolygons$Host)]

FullCentroids <-  FullPolygons %>% group_by(Host) %>% 
  summarise(LongMean = (max(long) + min(long))/2,
            LatMean = (max(lat) + min(lat))/2) %>% data.frame

rownames(FullCentroids) <- FullCentroids$Host

Panth1[,c("LongMean","LatMean")] <- FullCentroids[Panth1$Sp,c("LongMean","LatMean")]

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

#for(x in levels(Panth1$hOrder)){
#  hOrderList[[x]] <- induced_subgraph(AllSimGs[[1]], 
#                                      Panth1$hOrder == x)
#}

InDegree = lapply(hOrderList, function(a) lapply(a, function(b) degree(b)) %>% unlist)
InNames <- lapply(hOrderList[[1]], function(a) V(a)$name ) %>% unlist
InDegDF <- data.frame(InDegree)

InDegrees <- apply(InDegDF,1, function(a) mean(a, na.rm = T))

Panth1$InDegree <- InDegrees[Panth1$Sp]
Panth1$OutDegree <- with(Panth1, AllPredDegree-InDegree)

# Figs ####

BarGraph(Panth1, "hOrder", "AllPredDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  ggtitle("All Predicted Degrees") +#scale_fill_brewer(palette = AlberPalettes[2]) +
  ggsave("Figures/AllPredDegree.jpeg", units = "mm", height= 120, width = 200)

BarGraph(Panth1, "hOrder", "InDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  ggtitle("In-Degree") +#scale_fill_brewer(palette = AlberPalettes[2]) +
  ggsave("Figures/InDegree.jpeg", units = "mm", height= 120, width = 200)

BarGraph(Panth1, "hOrder", "OutDegree", order = T, text = "N") + 
  theme(legend.position = "none") + 
  ggtitle("Out-Degree") +#scale_fill_brewer(palette = AlberPalettes[2]) +
  ggsave("Figures/OutDegree.jpeg", units = "mm", height= 120, width = 200)


