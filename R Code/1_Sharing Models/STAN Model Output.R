
# STAN Model Output ####

library(rstan); library(reskew)

Prev <- function(x){
  
  length(x[x>0])/length(x)
  
}

logistic <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

traceplot(BinModel,
          pars = c("mu_alpha", 
                   "beta_d_cites_s", 
                   "beta_domestic", 
                   "beta_space",
                   "beta_phylo",
                   "beta_inter",
                   "sigma"))

p <- process_stanfit(BinModel, n.pars.to.trim = 3) # Takes a WHILE
#p <- process_stanfit(BinModel, n.pars.to.trim = 3) # Takes a WHILE

d = FinalHostMatrix

MCMCSol <- p$df#[rep(501:1000, 8)+rep(0:7*1000, each = 500),]

N = nrow(d)

d$Space_Phylo <- d$Space*d$Phylo2

d <- d %>% mutate(DCites = (log(hDiseaseZACites + 1)),#-mean(species.traits$d_cites))/sd(species.traits$d_cites), 
                  DCites.Sp2 = (log(hDiseaseZACites.Sp2 + 1)))#-mean(species.traits$d_cites))/sd(species.traits$d_cites))

# Simulating with specific random effects ####

XMatrix <- cbind(rep(1,N),
                 d[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

MZ1 <- model.matrix( ~ Sp - 1, data = FinalHostMatrix)
MZ2 <- model.matrix( ~ Sp2 - 1, data = FinalHostMatrix)

ZMatrixb <- MZ1 + MZ2 %>% as.matrix %>% as("dgCMatrix")
XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")

ZBetas2 <- colnames(p$df)[which(colnames(p$df)=="alpha_species[1]"):
                            which(colnames(p$df)=="alpha_species[649]")]

# Doing the simulating #####

PredList1 <- list()

RowsSampled <- sample(1:nrow(MCMCSol), 1000, replace = F)

for(x in 1:length(RowsSampled)){ # to do something non-specific
  
  if(x %% 10 == 0) print(x)
  
  XFX <- MCMCSol[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- MCMCSol[RowsSampled[x], ZBetas2] %>% unlist
  
  XPredictions <- XFX %*% t(XMatrix)
  
  ZPredictions <- ZFX %*% t(ZMatrixb)
  
  Predictions <- XPredictions + ZPredictions
  
  BinPred <- rbinom(n = N,
                    prob = logistic(Predictions@x),
                    size  = 1)
  
  PredList1[[x]] <- BinPred
  
}

PredDF1 <- data.frame(PredList1)

sapply(PredList1, Prev) %>% mean
sapply(PredList1, Prev) %>% qplot

d$PredVirus1 <- apply(PredDF1, 1, mean)

d$PredVirus1Q <- cut(d$PredVirus1,
                     breaks = c(-1:10/10),
                     labels = c(0:10/10))

# Simulating without random effects ####

XMatrix <- cbind(rep(1,N),
                 d[,c("Space","Phylo2","Space_Phylo")]) %>% as.matrix %>% as("dgCMatrix")

ZMatrix <- d[,c("DCites", "hDom","DCites.Sp2","hDom.Sp2")] %>% 
  mutate(hDom = ifelse(hDom == "wild", 0, 1), 
         hDom.Sp2 = ifelse(hDom.Sp2 == "wild", 0, 1)) %>%
  as.matrix %>% as("dgCMatrix")

XBetas <- c("mu_alpha","beta_space","beta_phylo","beta_inter")
ZBetas <- c("beta_d_cites_s","beta_domestic")

# Doing the simulating #####

PredList1b <- list()

RowsSampled <- sample(1:nrow(p$df), 1000, replace = F)

for(x in 1:length(RowsSampled)){ # to do something non-specific
  
  if(x %% 10 == 0) print(x)
  
  XFX <- p$df[RowsSampled[x], XBetas] %>% unlist
  
  ZFX <- p$df[RowsSampled[x], ZBetas] %>% unlist
  
  XPredictions <- XFX %*% t(XMatrix)
  
  ZPredictionsa <- ZFX %*% t(ZMatrix[,1:length(ZBetas)])
  ZPredictionsb <- ZFX %*% t(ZMatrix[,(length(ZBetas)+1):(length(ZBetas)*2)])
  
  ZPredictions2a <- rnorm(n = N, mean = ZPredictionsa@x, sd = p$df[RowsSampled[x], "sigma"])
  ZPredictions2b <- rnorm(n = N, mean = ZPredictionsb@x, sd = p$df[RowsSampled[x], "sigma"])
  #ZPredictions2c <- rnorm(n = N, mean = (ZPredictionsa@x + ZPredictionsb@x), sd = p$df[RowsSampled[x], "sigma"])
  
  Predictions <- XPredictions@x + ZPredictions2a + ZPredictions2b
  #Predictions2 <- XPredictions@x + ZPredictions2c
  
  BinPred <- rbinom(n = N,
                    prob = logistic(Predictions),
                    size  = 1)
  
  PredList1b[[x]] <- BinPred
  
}

PredDF1b <- data.frame(PredList1b)

sapply(PredList1b, Prev) %>% mean
sapply(PredList1b, Prev) %>% qplot

d$PredVirus1b <- apply(PredDF1b, 1, mean)

d$PredVirus1bQ <- cut(d$PredVirus1b,
                      breaks = c(-1:10/10),
                      labels = c(0:10/10))

# Simulating using the random effect ####

SimNets1 <- SimGraphs1 <- list()

for(i in 1:length(PredList1)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1[[i]] <- graph.adjacency(AssMat, mode = "undirected", diag = F)
  
}

Degdf1 <- sapply(SimGraphs1, function(a) degree(a))# %>% as.data.frame
Eigendf1 <- sapply(SimGraphs1, function(a) eigen_centrality(a)$vector) #%>% as.data.frame

PredDegrees1 <- apply(Degdf1, 1, mean)
PredEigen1 <- apply(Eigendf1, 1, mean)

Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]
Hosts$PredEigen1 <- PredEigen1[as.character(Hosts$Sp)]

# Simulating without the random effect ####

SimNets1b <- SimGraphs1b <- list()

for(i in 1:length(PredList1b)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- 0
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1b[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1b[[i]] <- graph.adjacency(AssMat, mode = "undirected")
  
}

Degdf1b <- sapply(SimGraphs1b, function(a) degree(a)) %>% as.data.frame
Eigendf1b <- sapply(SimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1b <- apply(Degdf1b, 1, mean)
PredEigen1b <- apply(Eigendf1b, 1, mean)

Hosts$PredDegree1b <- PredDegrees1b[as.character(Hosts$Sp)]
Hosts$PredEigen1b <- PredEigen1b[as.character(Hosts$Sp)]

ggplot(Hosts, aes(Degree, PredDegree1b)) + geom_point() + geom_smooth()

# Getting network-level stats

SimGraphList <- list(SimGraphs1, SimGraphs1b)# , SimGraphs2, SimGraphs3, SimGraphs3b)

Components <- lapply(SimGraphList, function(a) sapply(a, function(b) components(b)$no))
Degrees <- lapply(SimGraphList, function(a) sapply(a, function(b) mean(degree(b))))

Cluster1 = sapply(SimGraphs1, function(a) transitivity(a)) # all zero, don't bother
Cluster1b = sapply(SimGraphs1b, function(a) transitivity(a))

Betweenness1 = sapply(SimGraphs1, function(a) mean(betweenness(a)))
Betweenness1b = sapply(SimGraphs1b, function(a) mean(betweenness(a)))

Prev1 = apply(PredDF1, 2, Prev)
Prev1b = apply(PredDF1b, 2, Prev)

Closeness1 = sapply(SimGraphs1, function(a) mean(closeness(a)))
Closeness1b = sapply(SimGraphs1b, function(a) mean(closeness(a)))

