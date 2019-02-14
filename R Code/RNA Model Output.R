
# STAN Model Output ####

library(rstan); library(reskew); library(ggregplot); library(parallel)

if(file.exists("~/Albersnet/RNAModelOutput.Rdata")) load("~/Albersnet/RNAModelOutput.Rdata") else{
  
  RNABinModel <- readRDS("~/Albersnet/RNABinModel.rds")
  q <- process_stanfit(RNABinModel, n.pars.to.trim = 3) # Takes a WHILE
  save(q, file = "~/Albersnet/RNAModelOutput.Rdata")
}

MCMCSol <- q$df#[rep(501:1000, 8)+rep(0:7*1000, each = 500),]

e <- FinalHostMatrix %>% filter(!is.na(RNA))

e$Sp <- factor(as.character(e$Sp),
               levels = union(e$Sp, e$Sp2)
)

e$Sp2 <- factor(as.character(e$Sp2),
                levels = union(e$Sp, e$Sp2)
)

N = nrow(e)

UpperRNA <- which(upper.tri(matrix(NA, ncol = nlevels(e$Sp), nrow = nlevels(e$Sp2)), diag = T))

e$Space_Phylo <- scale(e$Space*e$Phylo2)

e$Space <- scale(e$Space)
e$Phylo2 <- scale(e$Phylo2)

e <- e %>% mutate(DCites = scale(log(hDiseaseZACites + 1)),#-mean(species.traits$d_cites))/sd(species.traits$d_cites), 
                  DCites.Sp2 = scale(log(hDiseaseZACites.Sp2 + 1)))#-mean(species.traits$d_cites))/sd(species.traits$d_cites))

# RNASimulating with specific random effects ####

XMatrix <- cbind(rep(1,N),
                 e[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

MZ1 <- model.matrix( ~ Sp - 1, data = e)
MZ2 <- model.matrix( ~ Sp2 - 1, data = e)

ZMatrixb <- MZ1 + MZ2 %>% as.matrix %>% as("dgCMatrix")
XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

ZBetas2 <- colnames(q$df)[which(colnames(q$df)=="alpha_species[1]"):
                            which(colnames(q$df)=="alpha_species[618]")]

# Doing the RNASimulating #####

RNAPredList1 <- list()

RowsSampled <- sample(1:nrow(MCMCSol), 1000, replace = F)

RNAPredList1 <- parallel::mclapply(1:length(RowsSampled), function(x){
  
  if(x %% 10 == 0) print(x)
  
  XFX <- MCMCSol[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- MCMCSol[RowsSampled[x], ZBetas2] %>% unlist
  
  XRNAPredictions <- XFX %*% t(XMatrix)
  
  ZRNAPredictions <- ZFX %*% t(ZMatrixb)
  
  RNAPredictions <- XRNAPredictions + ZRNAPredictions
  
  BinRNAPred <- rbinom(n = N,
                    prob = logistic(RNAPredictions@x),
                    size  = 1)
  
  BinRNAPred
  
}, mc.cores = 10)

RNAPredDF1 <- data.frame(RNAPredList1)
e$RNAPredVirus1 <- apply(RNAPredDF1, 1, mean)
e$RNAPredVirus1Q <- cut(e$RNAPredVirus1,
                     breaks = c(-1:10/10),
                     labels = c(0:10/10))

# RNASimulating without random effects ####

XMatrix <- cbind(rep(1,N),
                 e[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

ZMatrix <- e[,c("DCites", "hDom","DCites.Sp2","hDom.Sp2")] %>% 
  mutate(hDom = ifelse(hDom == "wild", 0, 1), 
         hDom.Sp2 = ifelse(hDom.Sp2 == "wild", 0, 1)) %>%
  as.matrix %>% as("dgCMatrix")

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

# Doing the RNASimulating #####

RNAPredList1b <- list()

RowsSampled <- sample(1:nrow(q$df), 1000, replace = F)

RNAPredList1b <- parallel::mclapply(1:length(RowsSampled), function(x){ # to do something non-specific
  
  XFX <- q$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- q$df[RowsSampled[x], ZBetas] %>% unlist
  
  XRNAPredictions <- XFX %*% t(XMatrix)
  
  ZRNAPredictionsa <- ZFX %*% t(ZMatrix[,1:length(ZBetas)])
  ZRNAPredictionsb <- ZFX %*% t(ZMatrix[,(length(ZBetas)+1):(length(ZBetas)*2)])
  
  ZRNAPredictions2a <- rnorm(n = N, mean = ZRNAPredictionsa@x, sd = q$df[RowsSampled[x], "sigma"])
  ZRNAPredictions2b <- rnorm(n = N, mean = ZRNAPredictionsb@x, sd = q$df[RowsSampled[x], "sigma"])
  #ZRNAPredictions2c <- rnorm(n = N, mean = (ZRNAPredictionsa@x + ZRNAPredictionsb@x), sd = q$df[RowsSampled[x], "sigma"])
  
  RNAPredictions <- XRNAPredictions@x + ZRNAPredictions2a + ZRNAPredictions2b
  #RNAPredictions2 <- XRNAPredictions@x + ZRNAPredictions2c
  
  BinRNAPred <- rbinom(n = N,
                    prob = logistic(RNAPredictions),
                    size  = 1)
  
  BinRNAPred
  
}, mc.cores = 10)

RNAPredDF1b <- data.frame(RNAPredList1b)
e$RNAPredVirus1b <- apply(RNAPredDF1b, 1, mean)
e$RNAPredVirus1bQ <- cut(e$RNAPredVirus1b,
                      breaks = c(-1:10/10),
                      labels = c(0:10/10))

# RNASimulating using the random effect ####

RNASimNets1 <- mclapply(1:length(RNAPredList1), function(i){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(e$Sp,e$Sp2)), 
                   ncol = length(union(e$Sp,e$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperRNA)] <- round(RNAPredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(e$Sp,e$Sp2),
                           union(e$Sp,e$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

RNASimGraphs1 <- mclapply(1:length(RNAPredList1), function(i){
  
  graph.adjacency(RNASimNets1[[i]], mode = "undirected", diag = F)
  
}, mc.cores = 10)

RNADegdf1 <- sapply(RNASimGraphs1, function(a) degree(a))# %>% as.data.frame
RNAPredDegrees1 <- apply(RNADegdf1, 1, mean)
Hosts$RNAPredDegree1 <- RNAPredDegrees1[as.character(Hosts$Sp)]

#Eigendf1 <- sapply(RNASimGraphs1, function(a) eigen_centrality(a)$vector) #%>% as.data.frame
#RNAPredEigen1 <- apply(Eigendf1, 1, mean)
#Hosts$RNAPredEigen1 <- RNAPredEigen1[as.character(Hosts$Sp)]

# RNASimulating without the random effect ####

RNASimNets1b <- mclapply(1:length(RNAPredList1b), function(i){
  
  
  AssMat <- matrix(NA, 
                   nrow = length(union(e$Sp,e$Sp2)), 
                   ncol = length(union(e$Sp,e$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperRNA)] <- round(RNAPredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(e$Sp,e$Sp2),
                           union(e$Sp,e$Sp2))
  
  as(AssMat, "dgCMatrix")
  
}, mc.cores = 10)

RNASimGraphs1b <- mclapply(1:length(RNAPredList1b), function(i){
  
  graph.adjacency(RNASimNets1b[[i]], mode = "undirected")
  
}, mc.cores = 10)

RNADegdf1b <- sapply(RNASimGraphs1b, function(a) degree(a)) %>% as.data.frame
RNAPredDegrees1b <- apply(RNADegdf1b, 1, mean)
Hosts$RNAPredDegree1b <- RNAPredDegrees1b[as.character(Hosts$Sp)]

#Eigendf1b <- sapply(RNASimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame
#RNAPredEigen1b <- apply(Eigendf1b, 1, mean)
#Hosts$RNAPredEigen1b <- RNAPredEigen1b[as.character(Hosts$Sp)]

# Getting network-level stats

RNASimGraphList <- list(RNASimGraphs1, RNASimGraphs1b)# , RNASimGraphs2, RNASimGraphs3, RNASimGraphs3b)

Components <- lapply(RNASimGraphList, function(a) sapply(a, function(b) components(b)$no))
ComponentSizes <- lapply(RNASimGraphList, function(a) sapply(a, function(b) components(b)$no))
lapply(RNASimGraphs1b[which(Components[[2]]==2)], function(a) which(components(a)$membership==2))

Degrees <- lapply(RNASimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(RNASimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1b = sapply(RNASimGraphs1b, function(a) transitivity(a))

#Betweenness1 = sapply(RNASimGraphs1, function(a) betweenness(a))
#Betweenness1b = sapply(RNASimGraphs1b, function(a) betweenness(a))

Prev1 = sapply(RNAPredList1, Prev)
Prev1b = sapply(RNAPredList1b, Prev)

#Closeness1 = sapply(RNASimGraphs1, function(a) closeness(a))
#Closeness1b = sapply(RNASimGraphs1b, function(a) closeness(a))

