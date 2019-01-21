
# Modelling Centrality Measures ####

library(INLA); library(tidyverse); library(ggregplot)

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

Resps <- c("Records", "Degree", "Eigenvector")

HostCentCovar = c(
  "hDom",
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  "PVRMass",
  #"hArtfclHbttUsrIUCN",
  #"hOrder",
  #"hMassGrams"
  "S.Greg1"
)

TestHosts <- Hosts %>% dplyr::select(Resps, HostCentCovar, "LongMean", "LatMean", "hOrder") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1), 
         Records = c(Records)) %>%
  
  slice(which(!NARows(Hosts[,c(Resps,HostCentCovar, "LongMean", "LatMean", "hOrder")])))


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
X <- as.data.frame(Xm[, -which(colnames(Xm) == "hDomdomestic")]) # Model Matrix

f1 = as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + ")))

f2 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)"))

f3 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

f4 =  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

FormulaList <- list(f1, f2, f3, f4)
ModelNames <- c("Base", "SPDE", "SpaceMat", "PDMat")
FamilyList <- c("poisson", "nbinomial", "beta")
CentralityList <- list()

for(r in 1:length(Resps)){ # Takes a while I bet
  
  CentralityList[[Resps[r]]] <- list()
  
  CentStack <- inla.stack(
    data = list(y = TestHosts[,Resp]),  
    A = list(1, 1, A3), # Vector of Multiplication factors              
    effects = list(
      Intercept = rep(1, N), # Leave
      X = X, # Leave
      w = w.index)) # Leave
  
  for(q in 1:length(ModelNames)){
    CentralityList[[Resps[r]]][[ModelNames[[q]]]] <-inla(
      FormulaList[[q]], 
      data = TestHosts, 
      family = "poisson",
      control.compute = list(dic = T))
  }
}


