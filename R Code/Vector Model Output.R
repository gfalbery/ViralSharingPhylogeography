
# STAN Model Output ####

library(rstan); library(reskew); library(ggregplot); library(parallel)

if(file.exists("~/Albersnet/VectorModelOutput.Rdata")) load("~/Albersnet/VectorModelOutput.Rdata") else{
  
  VectorBinModel <- readRDS("~/Albersnet/VectorBinModel.rds")
  s <- process_stanfit(VectorBinModel, n.pars.to.trim = 3) # Takes a WHILE
  save(s, file = "~/Albersnet/VectorModelOutput.Rdata")
}

g <- FinalHostMatrix %>% filter(!is.na(Vector))

g$Sp <- factor(as.character(g$Sp),
               levels = union(g$Sp, g$Sp2)
)

g$Sp2 <- factor(as.character(g$Sp2),
                levels = union(g$Sp, g$Sp2)
)

UpperVector <- which(upper.tri(matrix(NA, ncol = nlevels(g$Sp), nrow = nlevels(g$Sp2))))

MCMCSol <- s$df#[rep(501:1000, 8)+rep(0:7*1000, each = 500),]

N = nrow(g)

g$Space_Phylo <- scale(g$Space*g$Phylo2)
g$Space <- scale(g$Space)
g$Space_Phylo <- scale(g$Phylo2)

f <- f %>% mutate(DCites = (log(hDiseaseZACites + 1)),#-mean(species.traits$d_cites))/sd(species.traits$d_cites), 
                  DCites.Sp2 = (log(hDiseaseZACites.Sp2 + 1)))#-mean(species.traits$d_cites))/sd(species.traits$d_cites))

# VectorSimulating with specific random effects ####

XMatrix <- cbind(rep(1,N),
                 g[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

MZ1 <- model.matrix( ~ Sp - 1, data = f)
MZ2 <- model.matrix( ~ Sp2 - 1, data = f)

ZMatrixb <- MZ1 + MZ2 %>% as.matrix %>% as("dgCMatrix")
XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

ZBetas2 <- colnames(s$df)[which(colnames(s$df)=="alpha_species[1]"):
                            which(colnames(s$df)=="alpha_species[164]")]

# Doing the VectorSimulating #####

VectorPredList1 <- list()

RowsSampled <- sample(1:nrow(MCMCSol), 1000, replace = F)

VectorPredList1 <- parallel::mclapply(1:length(RowsSampled), function(x){
  
  if(x %% 10 == 0) print(x)
  
  XFX <- MCMCSol[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- MCMCSol[RowsSampled[x], ZBetas2] %>% unlist
  
  XVectorPredictions <- XFX %*% t(XMatrix)
  
  ZVectorPredictions <- ZFX %*% t(ZMatrixb)
  
  VectorPredictions <- XVectorPredictions + ZVectorPredictions
  
  BinVectorPred <- rbinom(n = N,
                       prob = logistic(VectorPredictions@x),
                       size  = 1)
  
  BinVectorPred
  
}, mc.cores = 10)

VectorPredDF1 <- data.frame(VectorPredList1)
g$VectorPredVirus1 <- apply(VectorPredDF1, 1, mean)
g$VectorPredVirus1Q <- cut(g$VectorPredVirus1,
                        breaks = c(-1:10/10),
                        labels = c(0:10/10))

# VectorSimulating without random effects ####

XMatrix <- cbind(rep(1,N),
                 g[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

ZMatrix <- g[,c("DCites", "hDom","DCites.Sp2","hDom.Sp2")] %>% 
  mutate(hDom = ifelse(hDom == "wild", 0, 1), 
         hDom.Sp2 = ifelse(hDom.Sp2 == "wild", 0, 1)) %>%
  as.matrix %>% as("dgCMatrix")

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

# Doing the VectorSimulating #####

VectorPredList1b <- list()

RowsSampled <- sample(1:nrow(s$df), 1000, replace = F)

VectorPredList1b <- parallel::mclapply(1:length(RowsSampled), function(x){ # to do something non-specific
  
  XFX <- s$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- s$df[RowsSampled[x], ZBetas] %>% unlist
  
  XVectorPredictions <- XFX %*% t(XMatrix)
  
  ZVectorPredictionsa <- ZFX %*% t(ZMatrix[,1:length(ZBetas)])
  ZVectorPredictionsb <- ZFX %*% t(ZMatrix[,(length(ZBetas)+1):(length(ZBetas)*2)])
  
  ZVectorPredictions2a <- rnorm(n = N, mean = ZVectorPredictionsa@x, sd = s$df[RowsSampled[x], "sigma"])
  ZVectorPredictions2b <- rnorm(n = N, mean = ZVectorPredictionsb@x, sd = s$df[RowsSampled[x], "sigma"])
  #ZVectorPredictions2c <- rnorm(n = N, mean = (ZVectorPredictionsa@x + ZVectorPredictionsb@x), sd = s$df[RowsSampled[x], "sigma"])
  
  VectorPredictions <- XVectorPredictions@x + ZVectorPredictions2a + ZVectorPredictions2b
  #VectorPredictions2 <- XVectorPredictions@x + ZVectorPredictions2c
  
  BinVectorPred <- rbinom(n = N,
                       prob = logistic(VectorPredictions),
                       size  = 1)
  
  BinVectorPred
  
}, mc.cores = 10)

VectorPredDF1b <- data.frame(VectorPredList1b)
g$VectorPredVirus1b <- apply(VectorPredDF1b, 1, mean)
g$VectorPredVirus1bQ <- cut(g$VectorPredVirus1b,
                         breaks = c(-1:10/10),
                         labels = c(0:10/10))

# VectorSimulating using the random effect ####

VectorSimNets1 <- mclapply(1:length(VectorPredList1), function(i){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(g$Sp,g$Sp2)), 
                   ncol = length(union(g$Sp,g$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperVector)] <- round(VectorPredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(g$Sp,g$Sp2),
                           union(g$Sp,g$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

VectorSimGraphs1 <- mclapply(1:length(VectorPredList1), function(i){
  
  graph.adjacency(VectorSimNets1[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

VectorDegdf1 <- sapply(VectorSimGraphs1, function(a) degree(a))# %>% as.data.frame
VectorPredDegrees1 <- apply(VectorDegdf1, 1, mean)
Hosts$VectorPredDegree1 <- VectorPredDegrees1[as.character(Hosts$Sp)]

#Eigendf1 <- sapply(VectorSimGraphs1, function(a) eigen_centrality(a)$vector) #%>% as.data.frame
#VectorPredEigen1 <- apply(Eigendf1, 1, mean)
#Hosts$VectorPredEigen1 <- VectorPredEigen1[as.character(Hosts$Sp)]

# VectorSimulating without the random effect ####

VectorSimNets1b <- mclapply(1:length(VectorPredList1b), function(i){
  
  
  AssMat <- matrix(NA, 
                   nrow = length(union(g$Sp,g$Sp2)), 
                   ncol = length(union(g$Sp,g$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperVector)] <- round(VectorPredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(g$Sp,g$Sp2),
                           union(g$Sp,g$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

VectorSimGraphs1b <- mclapply(1:length(VectorPredList1b), function(i){
  
  graph.adjacency(VectorSimNets1b[[i]], mode = "undirected")
  
}, mc.cores = 10)

VectorDegdf1b <- sapply(VectorSimGraphs1b, function(a) degree(a)) %>% as.data.frame
VectorPredDegrees1b <- apply(VectorDegdf1b, 1, mean)
Hosts$VectorPredDegree1b <- VectorPredDegrees1b[as.character(Hosts$Sp)]

Eigendf1b <- sapply(VectorSimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame
VectorPredEigen1b <- apply(Eigendf1b, 1, mean)
Hosts$VectorPredEigen1b <- VectorPredEigen1b[as.character(Hosts$Sp)]

# Getting network-level stats

VectorSimGraphList <- list(VectorSimGraphs1, VectorSimGraphs1b)# , VectorSimGraphs2, VectorSimGraphs3, VectorSimGraphs3b)

Components <- lapply(VectorSimGraphList, function(a) sapply(a, function(b) components(b)$no))
ComponentSizes <- lapply(VectorSimGraphList, function(a) sapply(a, function(b) components(b)$no))
lapply(VectorSimGraphs1b[which(Components[[2]]==2)], function(a) which(components(a)$membership==2))

Degrees <- lapply(VectorSimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(VectorSimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1b = sapply(VectorSimGraphs1b, function(a) transitivity(a))

#Betweenness1 = sapply(VectorSimGraphs1, function(a) betweenness(a))
#Betweenness1b = sapply(VectorSimGraphs1b, function(a) betweenness(a))

Prev1 = sapply(VectorPredList1, Prev)
Prev1b = sapply(VectorPredList1b, Prev)

#Closeness1 = sapply(VectorSimGraphs1, function(a) closeness(a))
#Closeness1b = sapply(VectorSimGraphs1b, function(a) closeness(a))

