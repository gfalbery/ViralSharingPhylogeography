
# Modelling Centrality Measures ####

# Making dataset ####

library(INLA); library(tidyverse); library(ggregplot)

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

Resps <- c("Records", "Degree", "Eigenvector", "hZoonosisCount")

HostCentCovar = c(
  "hDom",
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  #"hOrder",
  "S.Greg1"
)

TestHosts <- Hosts %>% dplyr::select(Resps, HostCentCovar, "Sp", "LongMean", "LatMean", "hOrder") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1), 
         Records = c(Records),
         Eigenvector = Eigenvector^0.2
  ) %>%
  
  slice(which(!NARows(Hosts[,c(Resps,HostCentCovar, "LongMean", "LatMean", "hOrder")])))  %>%
  dplyr::filter(!hOrder%in%c("SCANDENTIA","PERAMELEMORPHIA"))

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(HostCentCovar, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("hDomdomestic","hOrderCARNIVORA"))]) # Model Matrix

# Establishing distance matrices ####
SpaceContacts <- as(solve(RangeAdj1[FHN, FHN]),"dgCMatrix")
GRMatrix <- as(solve(tSTMatrix[FHN, FHN]),"dgCMatrix")

TestHosts$IndexSpace = unlist(sapply(TestHosts$Sp, function(a) which(FHN==a)))
TestHosts$IndexPhylo = unlist(sapply(TestHosts$Sp, function(a) which(FHN==a)))

# Establishing model formulae ####
f1 = as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + ")))

f2 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)"))

f3 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

f4 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

f5 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(w, model = spde)"))

f6 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(w, model = spde)"))

f7 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(IndexPhylo, model='generic0', Cmatrix = GRMatrix,
                       constr = TRUE,param = c(0.5, 0.5))"))

f8 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(IndexPhylo, model='generic0', Cmatrix = GRMatrix,
                       constr = TRUE,param = c(0.5, 0.5))",
                       "+ f(w, model = spde)"))

FormulaList <- list(f1, f2, f3, f4, f5, f6, f7)
ModelNames <- c("Base", "SPDE", "SpaceMat", "PDMat",
                "SpaceMat:w","PDMat:w","SpaceMat:PDMat")
FamilyList <- c("nbinomial", "nbinomial", "gaussian", "nbinomial")
CentralityList <- list()

for(r in 1:length(Resps)){ # Takes a while I bet
  
  print(Resps[r])
  
  CentralityList[[Resps[r]]] <- list()
  
  CentStack <- inla.stack(
    data = list(y = TestHosts[,Resps[r]]),  
    A = list(1, 1, 1, 1, A3), # Vector of Multiplication factors              
    effects = list(
      Intercept = rep(1, N), # Leave
      X = X, # Leave
      IndexSpace = TestHosts$IndexSpace,
      IndexPhylo = TestHosts$IndexPhylo,
      w = w.index)) # Leave
  
  for(q in 1:length(ModelNames)){
    print(ModelNames[q])
    CentralityList[[Resps[r]]][[ModelNames[[q]]]] <-
      inla(
        FormulaList[[q]], 
        family = FamilyList[r],
        data = inla.stack.data(CentStack),
        control.compute = list(dic = TRUE),
        control.predictor = list(A = inla.stack.A(CentStack))
      )
  }
  
  CentralityList[[Resps[r]]][[8]] <-
    inla(
      f8, 
      family = FamilyList[r],
      data = inla.stack.data(CentStack),
      control.compute = list(dic = TRUE),
      control.predictor = list(A = inla.stack.A(CentStack))
    )
  
}

# Plotting out ####

lapply(CentralityList, function(a) a$SPDE %>%
         ggField(WorldMesh)) %>% arrange_ggplot2(nrow = 3)

lapply(1:length(CentralityList), function(a) INLADICFig(CentralityList[[a]], ModelNames = ModelNames)) %>% 
  arrange_ggplot2(nrow = 3)

lapply(CentralityList, function(a) Efxplot(a, ModelNames = ModelNames)) %>% 
  arrange_ggplot2(ncol = 3)

lapply(1:length(CentralityList[[x]]), function(a){ 
  qplot(TestHosts[,Resps[x]],
        CentralityList[[x]][[a]]$summary.fitted.values$mean[1:dim(TestHosts)[1]]) + 
    ggtitle(ModelNames[a]) + labs(x = paste("Data", Resps[x]), y = paste("Fitted", Resps[x]))
}) %>%
  arrange_ggplot2

lapply(CentralityList, function(a) INLARep(a[[7]]))

# Plotting field ####

ggField(CentralityList[[3]][[2]], WorldMesh) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long/50000, lat/50000, group = group)) +
  geom_point(data = TestHosts[,c("LongMean", "LatMean")]/50000, aes(LongMean, LatMean), inherit.aes = F) + 
  scale_fill_brewer(palette = AlberPalettes[2])

hOrderCentralityList <- CentralityList

# Removing Order as a variable ####

HostCentCovar = c(
  "hDom",
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  #"PVRMass",
  #"hOrder",
  "S.Greg1"
)

TestHosts <- Hosts %>% dplyr::select(Resps, HostCentCovar, "Sp", "LongMean", "LatMean", "hOrder") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1), 
         Records = c(Records),
         Eigenvector = Eigenvector^0.2
  ) %>%
  
  slice(which(!NARows(Hosts[,c(Resps,HostCentCovar, "LongMean", "LatMean", "hOrder")])))  %>%
  dplyr::filter(!hOrder%in%c("SCANDENTIA","PERAMELEMORPHIA"))

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(HostCentCovar, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("hDomdomestic","hOrderCARNIVORA"))]) # Model Matrix

# Establishing model formulae ####
f1 = as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + ")))

f2 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)"))

f3 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

f4 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

f5 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(w, model = spde)"))

f6 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(w, model = spde)"))

f7 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       "+ f(IndexPhylo, model='generic0', Cmatrix = GRMatrix,
                       constr = TRUE,param = c(0.5, 0.5))"))

FormulaList <- list(f1, f2, f3, f4, f5, f6, f7)

# Establishing distance matrices ####
SpaceContacts <- as(solve(RangeAdj1[FHN, FHN]),"dgCMatrix")
GRMatrix <- as(solve(tSTMatrix[FHN, FHN]),"dgCMatrix")

TestHosts$IndexSpace = unlist(sapply(TestHosts$Sp, function(a) which(FHN==a)))
TestHosts$IndexPhylo = unlist(sapply(TestHosts$Sp, function(a) which(FHN==a)))

# Establishing model formulae ####
CentralityList <- list()

for(r in 1:length(Resps)){ # Takes a while I bet
  
  print(Resps[r])
  
  CentralityList[[Resps[r]]] <- list()
  
  CentStack <- inla.stack(
    data = list(y = TestHosts[,Resps[r]]),  
    A = list(1, 1, 1, 1, A3), # Vector of Multiplication factors              
    effects = list(
      Intercept = rep(1, N), # Leave
      X = X, # Leave
      IndexSpace = TestHosts$IndexSpace,
      IndexPhylo = TestHosts$IndexPhylo,
      w = w.index)) # Leave
  
  for(q in 1:length(ModelNames)){
    print(ModelNames[q])
    CentralityList[[Resps[r]]][[ModelNames[[q]]]] <-
      inla(
        FormulaList[[q]], 
        family = FamilyList[r],
        data = inla.stack.data(CentStack),
        control.compute = list(dic = TRUE),
        control.predictor = list(A = inla.stack.A(CentStack))
      )
  }
}

NohOrderCentralityList <- CentralityList

# Plotting out ####

lapply(CentralityList[1:4], function(a) a$SPDE %>%
         ggField(WorldMesh)) %>% arrange_ggplot2(nrow = 2)

lapply(1:length(CentralityList), function(a) INLADICFig(CentralityList[[a]])) %>% #, ModelNames = ModelNames)) %>% 
  arrange_ggplot2(nrow = 3)

lapply(CentralityList, function(a) Efxplot(a)) %>% #, ModelNames = ModelNames)) %>% 
  arrange_ggplot2(ncol = 3)

# Plotting field ####

ggField(CentralityList[[3]][[2]], WorldMesh) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long/50000, lat/50000, group = group)) +
  geom_point(data = TestHosts[,c("LongMean", "LatMean")]/50000, aes(LongMean, LatMean), inherit.aes = F) + 
  scale_fill_brewer(palette = AlberPalettes[2])
