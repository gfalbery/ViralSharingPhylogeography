
# STAN Model Output ####

library(rstan); library(reskew); library(ggregplot); library(parallel)

if(file.exists("~/Albersnet/NVectorModelOutput.Rdata")) load("~/Albersnet/NVectorModelOutput.Rdata") else{
  
  NVectorBinModel <- readRDS("~/Albersnet/NVectorBinModel.rds")
  t <- process_stanfit(NVectorBinModel, n.pars.to.trim = 3) # Takes a WHILE
  save(t, file = "~/Albersnet/NVectorModelOutput.Rdata")
}

h <- FinalHostMatrix %>% filter(!is.na(NVector))

h$Sp <- factor(as.character(h$Sp),
               levels = union(h$Sp, h$Sp2)
)

h$Sp2 <- factor(as.character(h$Sp2),
                levels = union(h$Sp, h$Sp2)
)

UpperNVector <- which(upper.tri(matrix(NA, ncol = nlevels(h$Sp), nrow = nlevels(h$Sp2))))

MCMCSol <- t$df#[rep(501:1000, 8)+rep(0:7*1000, each = 500),]

N = nrow(h)

h$Space_Phylo <- scale(h$Space*h$Phylo2)
h$Space <- scale(h$Space)
h$Space_Phylo <- scale(h$Phylo2)

f <- f %>% mutate(DCites = (log(hDiseaseZACites + 1)),#-mean(species.traitt$d_cites))/sd(species.traitt$d_cites), 
                  DCites.Sp2 = (log(hDiseaseZACites.Sp2 + 1)))#-mean(species.traitt$d_cites))/sd(species.traitt$d_cites))

# NVectorSimulating with specific random effects ####

XMatrix <- cbind(rep(1,N),
                 h[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

MZ1 <- model.matrix( ~ Sp - 1, data = f)
MZ2 <- model.matrix( ~ Sp2 - 1, data = f)

ZMatrixb <- MZ1 + MZ2 %>% as.matrix %>% as("dgCMatrix")
XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

ZBetas2 <- colnames(t$df)[which(colnames(t$df)=="alpha_species[1]"):
                            which(colnames(t$df)=="alpha_species[164]")]

# Doing the NVectorSimulating #####

NVectorPredList1 <- list()

RowsSampled <- sample(1:nrow(MCMCSol), 1000, replace = F)

NVectorPredList1 <- parallel::mclapply(1:length(RowsSampled), function(x){
  
  if(x %% 10 == 0) print(x)
  
  XFX <- MCMCSol[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- MCMCSol[RowsSampled[x], ZBetas2] %>% unlist
  
  XNVectorPredictions <- XFX %*% t(XMatrix)
  
  ZNVectorPredictions <- ZFX %*% t(ZMatrixb)
  
  NVectorPredictions <- XNVectorPredictions + ZNVectorPredictions
  
  BinNVectorPred <- rbinom(n = N,
                          prob = logistic(NVectorPredictions@x),
                          size  = 1)
  
  BinNVectorPred
  
}, mc.cores = 10)

NVectorPredDF1 <- data.frame(NVectorPredList1)
h$NVectorPredVirus1 <- apply(NVectorPredDF1, 1, mean)
h$NVectorPredVirus1Q <- cut(h$NVectorPredVirus1,
                           breaks = c(-1:10/10),
                           labels = c(0:10/10))

# NVectorSimulating without random effects ####

XMatrix <- cbind(rep(1,N),
                 h[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

ZMatrix <- h[,c("DCites", "hDom","DCites.Sp2","hDom.Sp2")] %>% 
  mutate(hDom = ifelse(hDom == "wild", 0, 1), 
         hDom.Sp2 = ifelse(hDom.Sp2 == "wild", 0, 1)) %>%
  as.matrix %>% as("dgCMatrix")

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

# Doing the NVectorSimulating #####

NVectorPredList1b <- list()

RowsSampled <- sample(1:nrow(t$df), 1000, replace = F)

NVectorPredList1b <- parallel::mclapply(1:length(RowsSampled), function(x){ # to do something non-specific
  
  XFX <- t$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- t$df[RowsSampled[x], ZBetas] %>% unlist
  
  XNVectorPredictions <- XFX %*% t(XMatrix)
  
  ZNVectorPredictionsa <- ZFX %*% t(ZMatrix[,1:length(ZBetas)])
  ZNVectorPredictionsb <- ZFX %*% t(ZMatrix[,(length(ZBetas)+1):(length(ZBetas)*2)])
  
  ZNVectorPredictions2a <- rnorm(n = N, mean = ZNVectorPredictionsa@x, sd = t$df[RowsSampled[x], "sigma"])
  ZNVectorPredictions2b <- rnorm(n = N, mean = ZNVectorPredictionsb@x, sd = t$df[RowsSampled[x], "sigma"])
  #ZNVectorPredictions2c <- rnorm(n = N, mean = (ZNVectorPredictionsa@x + ZNVectorPredictionsb@x), sd = t$df[RowsSampled[x], "sigma"])
  
  NVectorPredictions <- XNVectorPredictions@x + ZNVectorPredictions2a + ZNVectorPredictions2b
  #NVectorPredictions2 <- XNVectorPredictions@x + ZNVectorPredictions2c
  
  BinNVectorPred <- rbinom(n = N,
                          prob = logistic(NVectorPredictions),
                          size  = 1)
  
  BinNVectorPred
  
}, mc.cores = 10)

NVectorPredDF1b <- data.frame(NVectorPredList1b)
h$NVectorPredVirus1b <- apply(NVectorPredDF1b, 1, mean)
h$NVectorPredVirus1bQ <- cut(h$NVectorPredVirus1b,
                            breaks = c(-1:10/10),
                            labels = c(0:10/10))

# NVectorSimulating using the random effect ####

NVectorSimNets1 <- mclapply(1:length(NVectorPredList1), function(i){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(h$Sp,h$Sp2)), 
                   ncol = length(union(h$Sp,h$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperNVector)] <- round(NVectorPredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(h$Sp,h$Sp2),
                           union(h$Sp,h$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

NVectorSimGraphs1 <- mclapply(1:length(NVectorPredList1), function(i){
  
  graph.adjacency(NVectorSimNets1[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

NVectorDegdf1 <- sapply(NVectorSimGraphs1, function(a) degree(a))# %>% as.data.frame
NVectorPredDegrees1 <- apply(NVectorDegdf1, 1, mean)
Hostt$NVectorPredDegree1 <- NVectorPredDegrees1[as.character(Hosts$Sp)]

#Eigendf1 <- sapply(NVectorSimGraphs1, function(a) eigen_centrality(a)$NVector) #%>% as.data.frame
#NVectorPredEigen1 <- apply(Eigendf1, 1, mean)
#Hostt$NVectorPredEigen1 <- NVectorPredEigen1[as.character(Hostt$Sp)]

# NVectorSimulating without the random effect ####

NVectorSimNets1b <- mclapply(1:length(NVectorPredList1b), function(i){
  
  
  AssMat <- matrix(NA, 
                   nrow = length(union(h$Sp,h$Sp2)), 
                   ncol = length(union(h$Sp,h$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperNVector)] <- round(NVectorPredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(h$Sp,h$Sp2),
                           union(h$Sp,h$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

NVectorSimGraphs1b <- mclapply(1:length(NVectorPredList1b), function(i){
  
  graph.adjacency(NVectorSimNets1b[[i]], mode = "undirected")
  
}, mc.cores = 10)

NVectorDegdf1b <- sapply(NVectorSimGraphs1b, function(a) degree(a)) %>% as.data.frame
NVectorPredDegrees1b <- apply(NVectorDegdf1b, 1, mean)
Hostt$NVectorPredDegree1b <- NVectorPredDegrees1b[as.character(Hostt$Sp)]

Eigendf1b <- sapply(NVectorSimGraphs1b, function(a) eigen_centrality(a)$NVector) %>% as.data.frame
NVectorPredEigen1b <- apply(Eigendf1b, 1, mean)
Hostt$NVectorPredEigen1b <- NVectorPredEigen1b[as.character(Hostt$Sp)]

# Getting network-level stats

NVectorSimGraphList <- list(NVectorSimGraphs1, NVectorSimGraphs1b)# , NVectorSimGraphs2, NVectorSimGraphs3, NVectorSimGraphs3b)

Components <- lapply(NVectorSimGraphList, function(a) sapply(a, function(b) components(b)$no))
ComponentSizes <- lapply(NVectorSimGraphList, function(a) sapply(a, function(b) components(b)$no))
lapply(NVectorSimGraphs1b[which(Components[[2]]==2)], function(a) which(components(a)$membership==2))

Degrees <- lapply(NVectorSimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(NVectorSimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1b = sapply(NVectorSimGraphs1b, function(a) transitivity(a))

#Betweenness1 = sapply(NVectorSimGraphs1, function(a) betweenness(a))
#Betweenness1b = sapply(NVectorSimGraphs1b, function(a) betweenness(a))

Prev1 = sapply(NVectorPredList1, Prev)
Prev1b = sapply(NVectorPredList1b, Prev)

#Closeness1 = sapply(NVectorSimGraphs1, function(a) closeness(a))
#Closeness1b = sapply(NVectorSimGraphs1b, function(a) closeness(a))

