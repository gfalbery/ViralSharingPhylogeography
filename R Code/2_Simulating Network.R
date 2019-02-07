
# Generating the network of viruses/hosts

source("R Code/00_Master Code.R")

library(MCMCglmm); library(tidyverse)

invlogit <- function(a) exp(a)/(1 + exp(a))
logit <- function(a) log(a/(1-a))

load("Bin Model 1.Rdata")
load("Bin Model 2.Rdata")
load("~/Albersnet/Parallel_Binomials.Rdata")
load("~/Albersnet/Parallel_Binomialsb.Rdata")

# Checking Convergence ####

i = 1

mc <- BinModelList[1:10 + 10*(i-1)] %>% lapply(function(a) a$Sol[,1:14])
#mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$VCV)
mc <- do.call(mcmc.list, mc)
#par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
#gelman.plot(mc, auto.layout=F)
gelman.diag(mc)
par(mfrow=c(7,4), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

N = nrow(FinalHostMatrix)

# Going from the full model ####

PredList1 <- list()

ClusterMCMC <- BinModelList[1:10 + 10*(i-1)] %>% 
  lapply(function(a) as.data.frame(a$Sol)) %>% 
  bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- cbind(BinModelList[[i*10]]$X, BinModelList[[i*10]]$Z) %>% 
  as.matrix %>% as("dgCMatrix")

Columns <- list(1:ncol(BinModelList[[i*10]]$X),(ncol(BinModelList[[i*10]]$X)+1):ncol(XZMatrix))

for(x in 1:1000){
  if(x%%10==0) print(x)
  RowSampled <- RowsSampled[x]
  FXSample <- ClusterMCMC[RowSampled, unlist(Columns)]
  FXSample <- FXSample# + 0.04
  Output <- c(FXSample %*% t(XZMatrix))
  ProbVector <- Output[[1]]@x
  PZero <- rbinom(length(ProbVector), 1, logit(ProbVector))
  #PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector)*Prev(FinalHostMatrix$VirusBinary)/mean(invlogit(ProbVector)))
  PredList1[[x]] <- PZero
}

PredDF1 <- as.data.frame(PredList1)
names(PredDF1) <- paste("Rep",1:1000)
FinalHostMatrix$PredVirus1 <- apply(PredDF1,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus1Q <- cut(FinalHostMatrix$PredVirus1,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Simulating Networks ####

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
  
  SimGraphs1[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1 <- sapply(SimGraphs1, function(a) degree(a)) %>% as.data.frame
Eigendf1 <- sapply(SimGraphs1, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1 <- apply(Degdf1, 1, mean)
PredEigen1 <- apply(Eigendf1, 1, mean)

Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]
Hosts$PredEigen1 <- PredEigen1[as.character(Hosts$Sp)]

# Without random effects, same model ####

i = 1

PredList1b <- list()

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- BinModelList[[i*10]]$X %>% 
  as.matrix %>% as("dgCMatrix")

for(x in 1:1000){
  if(x%%10==0) print(x)
  RowSampled <- RowsSampled[x]
  FXSample <- ClusterMCMC[RowSampled, Columns[[1]]]
  Output <- c(FXSample %*% t(XZMatrix))
  ProbVector <- Output[[1]]@x
  PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector))
  #PZero <- rbinom(length(ProbVector), 1, invlogit(ProbVector)*Prev(FinalHostMatrix$VirusBinary)/mean(invlogit(ProbVector)))
  PredList1b[[x]] <- PZero
}

PredDF1b <- as.data.frame(PredList1b)
names(PredDF1b) <- paste("Rep",1:1000)
FinalHostMatrix$PredVirus1b <- apply(PredDF1b[,1:1000], 1, function(a) a %>% mean)

# Simulating the networks #####

SimNets1b <- SimGraphs1b <- list()

for(i in 1:length(PredList1b)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1b[[i]])# %>% as.matrix %>% as("dgCMatrix")
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))] #%>% as.matrix %>% as("dgCMatrix")
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1b[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs1b[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1b <- sapply(SimGraphs1b, function(a) degree(a)) %>% as.data.frame
Eigendf1b <- sapply(SimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1b <- apply(Degdf1b, 1, mean)
PredEigen1b <- apply(Eigendf1b, 1, mean)

Hosts$PredDegree1b <- PredDegrees1b[as.character(Hosts$Sp)]
Hosts$PredEigen1b <- PredEigen1b[as.character(Hosts$Sp)]

FinalHostMatrix$PredVirus1bQ <- cut(FinalHostMatrix$PredVirus1b,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Trying it without random effects ####

i = 2

PredList2 <- list()

ClusterMCMC <- BinModelList[1:10 + 10*(i-1)] %>% 
  lapply(function(a) as.data.frame(a$Sol)) %>% 
  bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- BinModelList[[1 + 10*(i-1)]]$X

for(x in 1:1000){
  if(x%%10==0) print(x)
  RowSampled <- RowsSampled[x]
  FXSample <- ClusterMCMC[RowSampled, Columns[[1]]]
  Output <- c(FXSample %*% t(XZMatrix))
  PZero <- rbinom(length(Output[[1]]@x), 1, invlogit(Output[[1]]@x))
  PredList2[[x]] <- PZero
}


PredDF2 <- as.data.frame(PredList2)
FinalHostMatrix$PredVirus2 <- apply(PredDF2,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus2Q <- cut(FinalHostMatrix$PredVirus2,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Simulating the network #####

SimNets2 <- SimGraphs2 <- list()

for(i in 1:length(PredList2)){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList2[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets2[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs2[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
  if(i%%10==0) print(i)
  
}

Degdf2 <- sapply(SimGraphs2, function(a) degree(a)) %>% as.data.frame
Eigendf2 <- sapply(SimGraphs2, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees2 <- apply(Degdf2, 1, mean)
PredEigen2 <- apply(Eigendf2, 1, mean)

Hosts$PredDegree2 <- PredDegrees2[as.character(Hosts$Sp)]
Hosts$PredEigen2 <- PredEigen2[as.character(Hosts$Sp)]

# Comparing differences, identifying which species are increased in the "true" ####

BeforeHPD <- cbind(apply(Degdf1[1:(nrow(Degdf1)/2),],1, function(a) HPDinterval(as.mcmc(a))[1]),
                   apply(Degdf1[1:(nrow(Degdf1)/2),],1, function(a) HPDinterval(as.mcmc(a))[2])) %>% as.data.frame

AfterHPD <- cbind(apply(Degdf1b[1:(nrow(Degdf1b)/2),],1, function(a) HPDinterval(as.mcmc(a))[1]),
                  apply(Degdf1b[1:(nrow(Degdf1b)/2),],1, function(a) HPDinterval(as.mcmc(a))[2])) %>% as.data.frame
BeforeHPD$When <- "Before"
AfterHPD$When <- "After"
BeforeHPD$Sp <- FHN
AfterHPD$Sp <- FHN
HPDComp <- rbind(BeforeHPD, AfterHPD)
HPDComp2 <- cbind(BeforeHPD, AfterHPD)
names(HPDComp2) <- paste(names(HPDComp2),rep(1:2,each = 4), sep = ".")

# Turning Up Citations with random effects ####

i = 1

PredList3 <- list()

ClusterMCMC <- BinModelList[1:10 + 10*(i-1)] %>% 
  lapply(function(a) as.data.frame(a$Sol)) %>% 
  bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- cbind(BinModelList[[i*10]]$X, BinModelList[[i*10]]$Z) %>% 
  as.matrix %>% as("dgCMatrix")

Columns <- list(1:ncol(BinModelList[[i*10]]$X),(ncol(BinModelList[[i*10]]$X)+1):ncol(XZMatrix))

XZMatrix[,"MinDCites"] <- max(XZMatrix[,"MinDCites"])
XZMatrix <- as(XZMatrix,"dgCMatrix")

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  FXSample <- ClusterMCMC[RowSampled, unlist(Columns)]
  Output <- c(FXSample %*% t(XZMatrix))
  
  PZero <- rbinom(length(Output[[1]]@x), 1, invlogit(Output[[1]]@x))
  
  PredList3[[x]] <- PZero
  
}

PredDF3 <- as.data.frame(PredList3)
FinalHostMatrix$PredVirus3 <- apply(PredDF3,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus3Q <- cut(FinalHostMatrix$PredVirus3,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Converting to graphs 

SimNets3 <- SimGraphs3 <- list()

for(i in 1:length(PredList3)){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList3[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets3[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs3[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
  if(i%%10==0) print(i)
  
}

Degdf3 <- sapply(SimGraphs3, function(a) degree(a)) %>% as.data.frame
Eigendf3 <- sapply(SimGraphs3, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees3 <- apply(Degdf3, 1, mean)
PredEigen3 <- apply(Eigendf3, 1, mean)

Hosts$PredDegree3 <- PredDegrees3[as.character(Hosts$Sp)]
Hosts$PredEigen3 <- PredEigen3[as.character(Hosts$Sp)]

# Turning Up Citations without random effects ####

PredList3b <- list()

ClusterMCMC <- BinModelList[1:10 + 10*(i-1)] %>% 
  lapply(function(a) as.data.frame(a$Sol)) %>% 
  bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- BinModelList[[10*i]]$X

XZMatrix[,"MinDCites"] <- max(XZMatrix[,"MinDCites"])
XZMatrix <- as(XZMatrix,"dgCMatrix")

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  FXSample <- ClusterMCMC[RowSampled, Columns[[1]]]
  Output <- c(FXSample %*% t(XZMatrix))
  
  PZero <- rbinom(length(Output[[1]]@x), 1, invlogit(Output[[1]]@x))
  
  PredList3b[[x]] <- PZero
  
}

PredDF3b <- as.data.frame(PredList3b)
FinalHostMatrix$PredVirus3b <- apply(PredDF3b,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus3bQ <- cut(FinalHostMatrix$PredVirus3b,
                                   breaks = c(-1:10/10),
                                   labels = c(0:10/10))

# Converting to graphs 

SimNets3b <- SimGraphs3b <- list()

for(i in 1:length(PredList3b)){
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList3b[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets3b[[i]] <- as(AssMat, "dgCMatrix")
  
  SimGraphs3b[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
  if(i%%10==0) print(i)
  
}

Degdf3b <- sapply(SimGraphs3b, function(a) degree(a)) %>% as.data.frame
Eigendf3b <- sapply(SimGraphs3b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees3b <- apply(Degdf3b, 1, mean)
PredEigen3b <- apply(Eigendf3b, 1, mean)

Hosts$PredDegree3b <- PredDegrees3b[as.character(Hosts$Sp)]
Hosts$PredEigen3b <- PredEigen3b[as.character(Hosts$Sp)]

# Looking at predicted changes in centrality ####

Hosts$DegreeChange <- with(Hosts, PredDegree3 - Degree)

BarGraph(Hosts,"hOrder","DegreeChange", text = "N")

# Comparing before and after ####

BeforeHPD <- cbind(apply(Degdf1[1:(nrow(Degdf1)/2),],1, function(a) HPDinterval(as.mcmc(a))[1]),
                   apply(Degdf1[1:(nrow(Degdf1)/2),],1, function(a) HPDinterval(as.mcmc(a))[2])) %>% as.data.frame

AfterHPD <- cbind(apply(Degdf1b[1:(nrow(Degdf3)/2),],1, function(a) HPDinterval(as.mcmc(a))[1]),
                  apply(Degdf1b[1:(nrow(Degdf3)/2),],1, function(a) HPDinterval(as.mcmc(a))[2])) %>% as.data.frame
BeforeHPD$When <- "Before"
AfterHPD$When <- "After"
BeforeHPD$Sp <- FHN
AfterHPD$Sp <- FHN
HPDComp <- rbind(BeforeHPD, AfterHPD)
HPDComp2 <- cbind(BeforeHPD, AfterHPD)
names(HPDComp2) <- paste(names(HPDComp2),rep(1:2,each = 4), sep = ".")

# Doing some checking ####

apply(FinalHostMatrix[,c("VirusBinary","PredVirus1","PredVirus1b","PredVirus2","PredVirus3")], 2, function(a) table(round(a, 2)==0))
apply(FinalHostMatrix[,c("VirusBinary","PredVirus1","PredVirus1b","PredVirus2","PredVirus3")], 2, function(a) mean(round(a, 2)))

FinalHostMatrix[,c("BinVirus1","BinVirus1b","BinVirus2","BinVirus3")] <- apply(FinalHostMatrix[,c("PredVirus1","PredVirus1b","PredVirus2","PredVirus3")],
                                                                               2, function(a) cut(a, breaks = c(-1,0:10/10), labels = c(0:10/10)))

# Trying an INLA model on the spatial distribution of degree changes ####

library(INLA)
TestHosts <- Hosts %>% dplyr::select("DegreeChange", "LongMean", "LatMean", "hOrder") %>%
  slice(which(!NARows(Hosts[,c("DegreeChange", "LongMean", "LatMean", "hOrder")])))  %>%
  dplyr::filter(!hOrder%in%c("SCANDENTIA","PERAMELEMORPHIA"))

N = nrow(TestHosts)

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Establishing model formulae ####
f1 = as.formula(paste("y ~ -1 + Intercept"))

f2 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        #paste(names(X), collapse = " + "),
                        "f(w, model = spde)"))

CentStack <- inla.stack(
  data = list(y = TestHosts$DegreeChange),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

IM1 <-
  inla(
    f1, 
    family = "gaussian",
    data = inla.stack.data(CentStack),
    control.compute = list(dic = TRUE),
    control.predictor = list(A = inla.stack.A(CentStack))
  )

IM2 <-
  inla(
    f2, 
    family = "gaussian",
    data = inla.stack.data(CentStack),
    control.compute = list(dic = TRUE),
    control.predictor = list(A = inla.stack.A(CentStack))
  )

INLADICFig(list(IM1, IM2))
ggField(IM2, WorldMesh)


