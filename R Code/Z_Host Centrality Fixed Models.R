
# Modelling Centrality Measures ####

library(INLA); library(tidyverse); library(ggregplot)

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

# Trying viral records as the response (similar to Olival et al., 2017) ####

Resp = "Records"

RecordsHostCentCovar = c(
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

TestHosts <- Hosts %>% dplyr::select(Resp, RecordsHostCentCovar) %>% #, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1), 
         #hMassGrams = log(hMassGrams),
         Records = c(Records)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,RecordsHostCentCovar, "LongMean", "LatMean")])))

IM1 <- INLAModelAdd(Resp, "1", RecordsHostCentCovar, Family = "nbinomial", Data = TestHosts)

RecordsKeptCovar <- RecordsHostCentCovar[-which(RecordsHostCentCovar%in%last(IM1$Removed))]

# Degree is the easiest to incorporate into network generation so:

Resp = "Degree"

DegreeHostCentCovar = c(
  "hAllZACites",
  #"hDiseaseZACites",
  #"hOrder",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  #"hMassGrams",
  "hDom"
)

TestHosts <- Hosts %>% dplyr::select(Resp, DegreeHostCentCovar, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange),
         #hMassGrams = log(hMassGrams), 
         S.Greg1 = sqrt(S.Greg1)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,DegreeHostCentCovar, "LongMean", "LatMean")])))

IM2 <- INLAModelAdd(Resp, "1", DegreeHostCentCovar, Family = "nbinomial", Data = TestHosts)

DegreeKeptCovar <- DegreeHostCentCovar[-which(DegreeHostCentCovar%in%last(IM2$Removed))]

# What if records is put as an explanatory variable? ####

Resp = "Degree"

DegreeHostCentCovar2 = c(
  "hAllZACites",
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  #"hMassGrams",
  "hDom",
  #"hOrder",
  "Records"
)

TestHosts <- Hosts %>% dplyr::select(Resp, DegreeHostCentCovar2, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         #hMassGrams = log(hMassGrams),
         Records = c(log(Records))) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,DegreeHostCentCovar2, "LongMean", "LatMean")])))

IM3 <- INLAModelAdd(Resp, "1", DegreeHostCentCovar2, Family = "nbinomial", Data = TestHosts)

DegreeKeptCovar2 <- DegreeHostCentCovar2[-which(DegreeHostCentCovar2%in%last(IM3$Removed))]

# Comparing the output to all these ####

Efxplot(lapply(list(IM1, IM2, IM3), function(a) a$FinalModel))

# Now trying Eigenvector j4j

Resp = "Eigenvector"

EigenHostCentCovar = c(
  "hAllZACites",
  #"hOrder"
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  #"hMassGrams",
  "hDom"
)

TestHosts <- Hosts %>% dplyr::select(Resp, EigenHostCentCovar, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange),
         #hMassGrams = log(hMassGrams), 
         S.Greg1 = sqrt(S.Greg1)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,EigenHostCentCovar, "LongMean", "LatMean")])))

IM4 <- INLAModelAdd(Resp, "1", EigenHostCentCovar, Family = "beta", Data = TestHosts)

EigenKeptCovar <- EigenHostCentCovar[-which(EigenHostCentCovar%in%last(IM4$Removed))]

Efxplot(lapply(list(IM1, IM2, IM3, IM4), function(a) a$FinalModel), 
        ModelNames = c("Records", "Raw Degree", "Record-Degree", "Eigen"))

# Eigenvector with Records included

Resp = "Eigenvector"

EigenHostCentCovar2 = c(
  "Records",
  "hAllZACites",
  #"hOrder"
  #"hDiseaseZACites",
  "GeogRange",
  "S.Greg1",
  "hArtfclHbttUsrIUCN",
  #"hMassGrams",
  "hDom"
)

TestHosts <- Hosts %>% dplyr::select(Resp, EigenHostCentCovar2, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange),
         #hMassGrams = log(hMassGrams), 
         S.Greg1 = sqrt(S.Greg1),
         Records = log(Records + 1)) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,EigenHostCentCovar2, "LongMean", "LatMean")])))

IM5 <- INLAModelAdd(Resp, "1", EigenHostCentCovar2, Family = "beta", Data = TestHosts)

EigenKeptCovar2 <- EigenHostCentCovar2[-which(EigenHostCentCovar2%in%last(IM5$Removed))]

# Adding Phylo or Space Matrices ####

tCytBMatrix <- 1 - (CytBMatrix - min(CytBMatrix))/max(CytBMatrix)

SpaceContacts <- as(solve(RangeAdj1[FHN, FHN]),"dgCMatrix")
GRMatrix <- as(solve(tCytBMatrix[FHN, FHN]),"dgCMatrix")

# NB phylo matrix has been inverted

FHosts <- Hosts[Hosts$Sp %in% FHN,]

FHosts$IndexSpace = unlist(sapply(FHosts$Sp, function(a) which(FHN==a)))
FHosts$IndexPhylo = unlist(sapply(FHosts$Sp, function(a) which(FHN==a)))

# Records ####

Resp = "Records"

TestHosts <- FHosts %>% dplyr::select(Resp, RecordsHostCentCovar, "LongMean", "LatMean", 
                                      "IndexSpace", "IndexPhylo") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         #hMassGrams = log(hMassGrams),
         Records = c(Records)) %>%
  
  slice(which(!NARows(FHosts[,c(Resp,RecordsHostCentCovar, "LongMean", "LatMean")])))

HostMMModelsRecords <- list()

f1 = as.formula(paste("Records ~ ", paste(RecordsKeptCovar, collapse = " + ")))

HostMMModelsRecords[[1]] <- inla(f1, 
                                 data = TestHosts, 
                                 family = "nbinomial",
                                 control.compute = list(dic = T))


f2 =  as.formula(paste("Records ~ ", paste(RecordsKeptCovar, collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsRecords[[2]] <- inla(f2, 
                                 data = TestHosts, 
                                 family = "nbinomial",
                                 control.compute = list(dic = T))


f3 =  as.formula(paste("Records ~ ", paste(RecordsKeptCovar, collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsRecords[[3]] <- inla(f3, 
                                 data = TestHosts, 
                                 family = "nbinomial",
                                 control.compute = list(dic = T))


# SPDE Model

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(RecordsKeptCovar, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) == "hDomdomestic")]) # Model Matrix

#f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
#f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid')"))

f4 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)"))

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    w = w.index)) # Leave

HostMMModelsRecords[[4]] <- inla(f4, # Base model (no random effects)
                             family = c("nbinomial"),
                             data = inla.stack.data(CentStack),
                             control.compute = list(dic = TRUE),
                             control.predictor = list(A = inla.stack.A(CentStack))
)

f4 = as.formula(paste("Records ~ ", paste(RecordsKeptCovar, collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))', 
                      " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsRecords[[4]] <- inla(f4, # Takes about an hour
                                 data = TestHosts, 
                                 family = "poisson",
                                 control.compute = list(dic = T))

sapply(HostMMModelsRecords, function(a) a$dic$dic)

# Degree ####

Resp = "Degree"

TestHosts <- FHosts %>% dplyr::select(Resp, DegreeKeptCovar, "LongMean", "LatMean", 
                                      "IndexSpace", "IndexPhylo") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         #hMassGrams = log(hMassGrams),
         Degree = c(Degree)) %>%
  
  slice(which(!NARows(FHosts[,c(Resp,DegreeHostCentCovar, "LongMean", "LatMean")])))

HostMMModelsDeg <- list()

f1 = as.formula(paste("Degree ~ ", paste(DegreeKeptCovar, collapse = " + ")))

HostMMModelsDeg[[1]] <- inla(f1, 
                             data = TestHosts, 
                             family = "poisson",
                             control.compute = list(dic = T))


f2 =  as.formula(paste("Degree ~ ", paste(DegreeKeptCovar, collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsDeg[[2]] <- inla(f2, 
                             data = TestHosts, 
                             family = "poisson",
                             control.compute = list(dic = T))


f3 =  as.formula(paste("Degree ~ ", paste(DegreeKeptCovar, collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsDeg[[3]] <- inla(f3, 
                             data = TestHosts, 
                             family = "poisson",
                             control.compute = list(dic = T))

# SPDE Model

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(DegreeKeptCovar, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) == "hDomdomestic")]) # Model Matrix

#f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
#f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid')"))

f4 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)"))

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    w = w.index)) # Leave

HostMMModelsDeg[[4]] <- inla(f4, # Base model (no random effects)
                             family = c("poisson"),
                             data = inla.stack.data(CentStack),
                             control.compute = list(dic = TRUE),
                             control.predictor = list(A = inla.stack.A(CentStack))
)

# Oh no get ready ####

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    IndexPhylo = TestHosts$IndexPhylo,
    w = w.index)) # Leave

HostMMModelsDeg[[5]] <- inla(f5, # Base model (no random effects)
                               family = c("beta"),
                               data = inla.stack.data(CentStack),
                               control.compute = list(dic = TRUE),
                               control.predictor = list(A = inla.stack.A(CentStack))
)

# Multi-Matrix Model

f4 = as.formula(paste("Degree ~ ", paste(DegreeKeptCovar, collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))', 
                      " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsDeg[[4]] <- inla(f4,
                             data = TestHosts, 
                             family = "poisson",
                             control.compute = list(dic = T))

sapply(HostMMModelsDeg, function(a) a$dic$dic)

# Eigenvector #####

Resp = "Eigenvector"

TestHosts <- FHosts %>% dplyr::select(Resp, EigenKeptCovar, "LongMean", "LatMean", 
                                      "IndexSpace", "IndexPhylo") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange), 
         S.Greg1 = sqrt(S.Greg1),
         #hMassGrams = log(hMassGrams),
         Eigenvector = c(Eigenvector)) %>%
  
  slice(which(!NARows(FHosts[,c(Resp,DegreeHostCentCovar, "LongMean", "LatMean")])))

HostMMModelsEigen <- list()

f1 = as.formula(paste("Eigenvector ~ ", paste(EigenKeptCovar, collapse = " + ")))

HostMMModelsEigen[[1]] <- inla(f1, 
                               data = TestHosts, 
                               family = "beta",
                               control.compute = list(dic = T))


f2 =  as.formula(paste("Eigenvector ~ ", paste(EigenKeptCovar, collapse = " + "), " + ", 'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsEigen[[2]] <- inla(f2, 
                               data = TestHosts, 
                               family = "beta",
                               control.compute = list(dic = T))


f3 =  as.formula(paste("Eigenvector ~ ", paste(EigenKeptCovar, collapse = " + "), " + ", 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

HostMMModelsEigen[[3]] <- inla(f3, 
                               data = TestHosts, 
                               family = "beta",
                               control.compute = list(dic = T))

# SPDE Model

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000

HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)

WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)

A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix

spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(EigenKeptCovar, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("hDomdomestic"))]) # Model Matrix

#f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
#f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid')"))

f4 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)"))

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    w = w.index)) # Leave

HostMMModelsEigen[[4]] <- inla(f4, # Base model (no random effects)
                               family = c("beta"),
                               data = inla.stack.data(CentStack),
                               control.compute = list(dic = TRUE),
                               control.predictor = list(A = inla.stack.A(CentStack))
)

sapply(HostMMModelsEigen, function(a) a$dic$dic)

ggField(HostMMModelsEigen[[4]], WorldMesh) + geom_point(data = TestHosts, inherit.aes = F, aes(x = LongMean, y = LatMean))

# Oh No Get Ready for the phylo matrix + space too wahey ####

f5 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)",
                        ' + f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'))

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    IndexPhylo = TestHosts$IndexPhylo,
    w = w.index)) # Leave

HostMMModelsEigen[[5]] <- inla(f5, # Base model (no random effects)
                               family = c("beta"),
                               data = inla.stack.data(CentStack),
                               control.compute = list(dic = TRUE),
                               control.predictor = list(A = inla.stack.A(CentStack))
)

# Oh no it did not help ####

sapply(HostMMModelsEigen, function(a) a$dic$dic)
Efxplot(HostMMModelsEigen)

# Trying model addition to model with phylo ####

Expl = 'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix, constr=TRUE, param = c(0.5, 0.5))'

Expl = 'f(IndexSpace, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))'

Resps = c("Records", "Degree", "Eigenvector")

TestHosts <- FHosts %>% dplyr::select(Resps, DegreeKeptCovar, "LongMean", "LatMean", "IndexPhylo") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1),
         #hDiseaseZACites = log(hDiseaseZACites + 1),
         GeogRange = kader:::cuberoot(GeogRange),
         #hMassGrams = log(hMassGrams), 
         S.Greg1 = sqrt(S.Greg1)) %>%
  
  slice(which(!NARows(Hosts[,c(Resps,DegreeKeptCovar, "LongMean", "LatMean")])))

FullAdd1 <- INLAModelAdd(Resps[1], Expl, RecordsKeptCovar, Family = "nbinomial", Data = TestHosts)
FullAdd2 <- INLAModelAdd(Resps[2], Expl, DegreeKeptCovar, Family = "poisson", Data = TestHosts)
FullAdd3 <- INLAModelAdd(Resps[3], Expl, EigenKeptCovar, Family = "beta", Data = TestHosts)


