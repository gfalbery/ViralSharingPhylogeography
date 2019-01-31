
# Generating the network of viruses/hosts

library(MCMCglmm); library(tidyverse)

load("ZI_runs.Rdata")

CountColumns <- list(1:7*2-1, 15:662)
ZIColumns <- list(1:7*2, 663:1310)

i = 1

N = nrow(FinalHostMatrix)

PredList1 <- list()

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

logit <- function(a) exp(a)/(1 + exp(a))

XZMatrix <- cbind(ZI_runs[[(i-1)*10+1]]$X, ZI_runs[[(i-1)*10+1]]$Z)

CountXZMatrix <- XZMatrix[1:N,unlist(CountColumns)] #%>% as.matrix 
ZIXZMatrix <- XZMatrix[(N+1):(2*N),unlist(ZIColumns)] #%>% as.matrix

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, unlist(CountColumns)]
  ZIFXSample <- ClusterMCMC[RowSampled, unlist(ZIColumns)]
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  #PZero <- logit(ZIOutput[[1]]@x)
  #PCount <- exp(CountOutput[[1]]@x)*(1-PZero)
  
  PZero <- rbinom(length(ZIOutput[[1]]@x), 1, logit(ZIOutput[[1]]@x))
  PCount <- rpois(length(ZIOutput[[1]]@x),exp(CountOutput[[1]]@x))*(1-PZero)
  
  PCount <- ifelse(PCount>0,1,0)
  
  PredList1[[x]] <- PCount
  
}

PredDF1 <- as.data.frame(PredList1)
names(PredDF1) <- paste("Rep",1:1000)
PredDF1$Actual <- FinalHostMatrix$VirusBinary
lapply(1:length(HPD), function(b) FinalHostMatrix)

MeanPredictions <- apply(PredDF1,1, function(a) a %>% mean)
ModePredictions <- apply(PredDF1,1, function(a) a %>% median)

FinalHostMatrix$PredVirus1 <- ModePredictions

# Without random effects, same model ####

PredList1b <- list()

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

CountXZMatrix <- XZMatrix[1:N,CountColumns[[1]]] #%>% as.matrix 
ZIXZMatrix <- XZMatrix[(N+1):(2*N),ZIColumns[[1]]] #%>% as.matrix

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, CountColumns[[1]]]
  ZIFXSample <- ClusterMCMC[RowSampled, ZIColumns[[1]]]
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  #PZero <- logit(ZIOutput[[1]]@x)
  #PCount <- exp(CountOutput[[1]]@x)*(1-PZero)
  
  PZero <- rbinom(length(ZIOutput[[1]]@x), 1, logit(ZIOutput[[1]]@x))
  PCount <- rpois(length(ZIOutput[[1]]@x),exp(CountOutput[[1]]@x))*(1-PZero)
  
  PCount <- ifelse(PCount>0,1,0)
  
  PredList1b[[x]] <- PCount
  
}

PredDF1b <- as.data.frame(PredList1b)
MeanPredictions <- apply(PredDF1b,1, function(a) a %>% mean)
FinalHostMatrix$PredVirus1b <- MeanPredictions

# Trying it without random effects ####

i = 2

PredList2 <- list()

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

XZMatrix <- cbind(ZI_runs[[(i-1)*10+1]]$X, ZI_runs[[(i-1)*10+1]]$Z)

CountXZMatrix <- XZMatrix[1:N,CountColumns[[1]]] #%>% as.matrix
ZIXZMatrix <- XZMatrix[(N+1):(2*N),ZIColumns[[1]]] #%>% as.matrix

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, CountColumns[[1]]]
  ZIFXSample <- ClusterMCMC[RowSampled, ZIColumns[[1]]]
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput[[1]]@x)
  PCount <- exp(CountOutput[[1]]@x)*(1-PZero)
  PCount <- ifelse(PCount>0,1,0)
  
  PredList2[[x]] <- PCount
  
}

PredDF2 <- as.data.frame(PredList2)
ModePredictions <- apply(PredDF2,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus2 <- ModePredictions

GGally::ggpairs(FinalHostMatrix[,c("Virus","PredVirus1","PredVirus1b","PredVirus2")],
                lower = list(continuous = "smooth"))

# Simulating the networks #####

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

Degdflong1 <- reshape2::melt(t(Degdf1)) %>% rename(Sp = Var2, Degree = value)

a = 850

ggplot(Degdflong1[1:a,], aes(Sp, Degree)) + geom_violin(aes(colour = Sp)) +
  geom_point(data = Hosts[Hosts$Sp%in%Degdflong1[1:a,"Sp"],], aes(Sp, Degree)) + 
  facet_wrap(~Sp, scales = "free")

PredDegrees1 <- apply(Degdf1, 1, mean)

Hosts$PredDegree1 <- PredDegrees1[as.character(Hosts$Sp)]

# Simulating the networks #####

SimNets1b <- SimGraphs1b <- list()

for(i in 1:length(PredList1b)){
  
  if(i%%10==0) print(i)
  
  AssMat <- matrix(NA, 
                   nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                   ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
  
  AssMat[-which(1:length(AssMat)%in%UpperHosts)] <- round(PredList1b[[i]])
  AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
  diag(AssMat) <- apply(AssMat,1,function(a) length(a[!is.na(a)&a>0]))
  dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                           union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
  
  SimNets1b[[i]] <- AssMat
  
  SimGraphs1b[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
}

Degdf1b <- sapply(SimGraphs1b, function(a) degree(a)) %>% as.data.frame
Eigendf1b <- sapply(SimGraphs1b, function(a) eigen_centrality(a)$vector) %>% as.data.frame

PredDegrees1b <- apply(Degdf1b, 1, mean)

Hosts$PredDegree1b <- PredDegrees1b[as.character(Hosts$Sp)]

# Trying sans random effect ####

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
  
  SimNets2[[i]] <- AssMat
  
  SimGraphs2[[i]] <- graph.incidence(AssMat, weighted = TRUE)
  
  if(i%%10==0) print(i)
  
}

Degdf2 <- sapply(SimGraphs2, function(a) degree(a)) %>% as.data.frame
Eigendf2 <- sapply(SimGraphs2[1:729], function(a) eigen_centrality(a)$vector) %>% as.data.frame

Degdflong2 <- reshape2::melt(t(Degdf2)) %>% rename(Sp = Var2, Degree = value)

PredDegrees2 <- apply(Degdf2, 1, mean)

Hosts$PredDegree2 <- PredDegrees2[as.character(Hosts$Sp)]

ggplot(Hosts, aes(PredDegree1b, PredDegree2)) + geom_point() + geom_smooth()

GGally::ggpairs(Hosts[,c("Degree","PredDegree1","PredDegree1b","PredDegree2")],
                lower = list(continuous = "smooth"))

apply(Hosts[,c("Degree","PredDegree1","PredDegree1b","PredDegree2")],2,function(a) table(a>0))

PredEigen1 <- apply(Eigendf1,1,mean)
PredEigen1b <- apply(Eigendf1b,1,mean)
PredEigen2 <- apply(Eigendf2,1,mean)

Hosts[,c("Eigen1","Eigen1b","Eigen2")] <- cbind(PredEigen1[as.character(Hosts$Sp)],
                                                PredEigen1b[as.character(Hosts$Sp)],
                                                PredEigen2[as.character(Hosts$Sp)])

GGally::ggpairs(Hosts[,c("Eigenvector","Eigen1","Eigen1b","Eigen2")],
                lower = list(continuous = "smooth"))


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

# Turning Up Citations ####

i = 1

PredList3 <- list()

XZMatrix <- cbind(ZI_runs[[(i-1)*10+1]]$X, ZI_runs[[(i-1)*10+1]]$Z)

CountXZMatrix <- XZMatrix[1:N,unlist(CountColumns)] %>% as.matrix
ZIXZMatrix <- XZMatrix[(N+1):(2*N),unlist(ZIColumns)] %>% as.matrix

CountXZMatrix[,"traitVirus:MinDCites"] <- max(CountXZMatrix[,"traitVirus:MinDCites"])
ZIXZMatrix[,"traitzi_Virus:MinDCites"] <- max(ZIXZMatrix[,"traitzi_Virus:MinDCites"])

CountXZMatrix <- as(CountXZMatrix,"dgCMatrix")
ZIXZMatrix <- as(ZIXZMatrix,"dgCMatrix")

N = nrow(FinalHostMatrix)

ClusterMCMC <- ZI_runs[1:10 + (i-1)*10] %>% lapply(function(a) as.data.frame(as.matrix(a$Sol))) %>% bind_rows %>% as.matrix

RowsSampled <- sample(1:nrow(ClusterMCMC), 1000, replace = F)

for(x in 1:1000){
  
  if(x%%10==0) print(x)
  
  RowSampled <- RowsSampled[x]
  
  CountFXSample <- ClusterMCMC[RowSampled, unlist(CountColumns)]
  ZIFXSample <- ClusterMCMC[RowSampled, unlist(ZIColumns)]
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- rbinom(length(ZIOutput[[1]]@x), 1, logit(ZIOutput[[1]]@x))
  PCount <- rpois(length(ZIOutput[[1]]@x),exp(CountOutput[[1]]@x))*(1-PZero)
  
  PCount <- ifelse(PCount>0,1,0)
  
  PredList3[[x]] <- PCount
  
}

PredDF3 <- as.data.frame(PredList3)
ModePredictions <- apply(PredDF3,1, function(a) a %>% mean)

FinalHostMatrix$PredVirus3 <- ModePredictions

FinalHostMatrix %>% ggplot(aes(VirusBinary, PredVirus3)) + geom_point() + geom_smooth()

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
  
  SimGraphs3[[i]] <- graph.incidence(AssMat)
  
  if(i%%10==0) print(i)
  
}

Degdf3 <- sapply(SimGraphs3, function(a) degree(a)) %>% as.data.frame
Eigendf3 <- sapply(SimGraphs3, function(a) eigen_centrality(a)$vector) %>% as.data.frame

Degdflong3 <- reshape2::melt(t(Degdf3)) %>% rename(Sp = Var2, Degree = value)

PredDegrees3 <- apply(Degdf3, 1, mean)

Hosts$PredDegree3 <- PredDegrees3[as.character(Hosts$Sp)]

ggplot(Hosts,aes(Degree, PredDegree3)) + geom_point() + geom_smooth()

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

apply(FinalHostMatrix[,c("Virus","PredVirus1","PredVirus1b","PredVirus2","PredVirus3")], 2, function(a) table(round(a, 2)==0))
apply(FinalHostMatrix[,c("Virus","PredVirus1","PredVirus1b","PredVirus2","PredVirus3")], 2, function(a) mean(round(a, 2)))

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


