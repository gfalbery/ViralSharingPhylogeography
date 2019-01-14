
# Fitting Animal Models using matrices ####

# Making Host Data ####

SpaceContacts <- as(solve(RangeAdj1[FHN, FHN]),"dgCMatrix")
GRMatrix <- as(solve(1-CytBMatrix[FHN, FHN]),"dgCMatrix")

# NB phylo matrix has been inverted

FHosts$IndexSpace = sapply(FHosts$Sp, function(a) which(FHN==a))
FHosts$IndexPhylo = sapply(FHosts$Sp, function(a) which(FHN==a))

TestHosts <- FHosts

TestHosts$hEigenvector <- kader:::cuberoot(TestHosts$Eigenvector)

HostMMModels <- list()

f1 = hEigenvector ~ 1

HostMMModels[[1]] <- inla(f1, 
                          data = TestHosts, 
                          family = "gaussian",
                          control.compute = list(dic = T))


f2 =  hEigenvector ~ f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))

HostMMModels[[2]] <- inla(f2, 
                          data = TestHosts, 
                          family = "gaussian",
                          control.compute = list(dic = T))


f3 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5))

HostMMModels[[3]] <- inla(f3, 
                          data = TestHosts, 
                          family = "gaussian",
                          control.compute = list(dic = T))

f4 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

HostMMModels[[4]] <- inla(f4, # Takes about an hour ####
                          data = TestHosts, 
                          family = "gaussian",
                          control.compute = list(dic = T))

sapply(HostMMModels, function(a) a$dic$dic)

save(HostMMModels, file = "Model Files/HostMMModels.Rdata")

# Fitting space sharing and phylogenies both improve centrality fit 

# Extracting variance readings ####

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  HostMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace = 1/inla.qmarginal(c(0.025, 0.975), 
                                 HostMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo = inla.emarginal(function(x) 1/x,
                                  HostMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo = 1/inla.qmarginal(c(0.025, 0.975), 
                                 HostMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev

# From joint models ####

sigma.IndexSpace2 = inla.emarginal(function(x) 1/x,
                                   HostMMModels[[4]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace2 = 1/inla.qmarginal(c(0.025, 0.975), 
                                  HostMMModels[[4]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo2 = inla.emarginal(function(x) 1/x,
                                   HostMMModels[[4]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo2 = 1/inla.qmarginal(c(0.025, 0.975), 
                                  HostMMModels[[4]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev

# Making Virus dataset ####

VirusSpaceContacts <- as(solve(VirusRangeAdj1[FVN2, FVN2] + matrix(runif(length(FVN2)^2, min=0,max=0.01), ncol = length(FVN2))),"dgCMatrix")
VirusGRMatrix <- as(solve(1-VirusHostPD[FVN2, FVN2] + matrix(runif(length(FVN2)^2, min=0,max=0.01), ncol = length(FVN2))),"dgCMatrix")

# NB phylo matrix has been inverted ####

FViruses <- FViruses[FViruses$Sp%in%FVN2,]

FViruses$IndexSpace = sapply(FViruses$Sp, function(a) which(FVN2 == a))
FViruses$IndexPhylo = sapply(FViruses$Sp, function(a) which(FVN2 == a))

TestViruses <- FViruses

TestViruses$vEigenvector <- kader:::cuberoot(TestViruses$Eigenvector)

VirusMMModels <- list()

f1 = vEigenvector ~ 1

VirusMMModels[[1]] <- inla(f1, 
                           data = TestViruses, 
                           family = "gaussian",
                           control.compute = list(dic = T))

f2 =  vEigenvector ~ f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
                       constr=TRUE, param = c(0.5, 0.5))

VirusMMModels[[2]] <- inla(f2, 
                           data = TestViruses, 
                           family = "gaussian",
                           control.compute = list(dic = T))

f3 = vEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = VirusGRMatrix,
                      constr = TRUE, param = c(0.5, 0.5))

VirusMMModels[[3]] <- inla(f3, 
                           data = TestViruses, 
                           family = "gaussian",
                           control.compute = list(dic = T))

f4 = vEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = VirusGRMatrix,
                      constr = TRUE, param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

VirusMMModels[[4]] <- inla(f4, # This takes a long time to run 
                           data = TestViruses, 
                           family = "gaussian",
                           control.compute = list(dic = T))

sapply(VirusMMModels, function(a) a$dic$dic)

# Fitting matrices do not improve centrality model fit #### 

# Extracting variance readings ####

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  VirusMMModels[[2]]$marginals.hyperpar$"Precision for IndexSpace")

e.IndexA = inla.expectation(function(x) x, sigma.IndexSpace)
ci.IndexA = inla.qmarginal(c(0.025, 0.975), sigma.IndexSpace)

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  VirusMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace = 1/inla.qmarginal(c(0.025, 0.975), 
                                 VirusMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo = inla.emarginal(function(x) 1/x,
                                  VirusMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo = 1/inla.qmarginal(c(0.025, 0.975), 
                                 VirusMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev

# Do space and phylogeny affect zoonotic capacity? ####

ZooMMModels <- list()

f1 = Human ~ 1

ZooMMModels[[1]] <- inla(f1, 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

f2 =  Human ~ f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
                constr=TRUE, param = c(0.5, 0.5))

ZooMMModels[[2]] <- inla(f2, 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

f3 = Human ~ f(IndexPhylo, model = "generic0", Cmatrix = VirusGRMatrix,
               constr = TRUE, param = c(0.5, 0.5))

ZooMMModels[[3]] <- inla(f3, 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

f4 = Human ~ f(IndexPhylo, model = "generic0", Cmatrix = VirusGRMatrix,
               constr = TRUE, param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

ZooMMModels[[4]] <- inla(f4, # This takes a long time to run 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

sapply(ZooMMModels, function(a) a$dic$dic)

# Nothing improves zoonotic 0/1 models 

# Extracting variance readings ####

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  ZooMMModels[[2]]$marginals.hyperpar$"Precision for IndexSpace")

e.IndexA = inla.expectation(function(x) x, sigma.IndexSpace)
ci.IndexA = inla.qmarginal(c(0.025, 0.975), sigma.IndexSpace)

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  ZooMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace = 1/inla.qmarginal(c(0.025, 0.975), 
                                 ZooMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo = inla.emarginal(function(x) 1/x,
                                  ZooMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo = 1/inla.qmarginal(c(0.025, 0.975), 
                                 ZooMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev

# Do space and phylogeny affect domestic infection capacity? ####

DomMMModels <- list()

f1 = Domestic ~ 1

DomMMModels[[1]] <- inla(f1, 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

f2 =  Domestic ~ f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
                constr=TRUE, param = c(0.5, 0.5))

DomMMModels[[2]] <- inla(f2, 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

f3 = Domestic ~ f(IndexPhylo, model = "generic0", Cmatrix = VirusGRMatrix,
               constr = TRUE, param = c(0.5, 0.5))

DomMMModels[[3]] <- inla(f3, 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

f4 = Domestic ~ f(IndexPhylo, model = "generic0", Cmatrix = VirusGRMatrix,
               constr = TRUE, param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

DomMMModels[[4]] <- inla(f4, # This takes a long time to run 
                         data = TestViruses, 
                         family = "binomial",
                         control.compute = list(dic = T))

sapply(DomMMModels, function(a) a$dic$dic)

# Fitting space sharing and phylogenies both improve centrality fit 

# Extracting variance readings ####

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  DomMMModels[[2]]$marginals.hyperpar$"Precision for IndexSpace")

e.IndexA = inla.expectation(function(x) x, sigma.IndexSpace)
ci.IndexA = inla.qmarginal(c(0.025, 0.975), sigma.IndexSpace)

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  DomMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace = 1/inla.qmarginal(c(0.025, 0.975), 
                                 DomMMModels[[2]]["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo = inla.emarginal(function(x) 1/x,
                                  DomMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo = 1/inla.qmarginal(c(0.025, 0.975), 
                                 DomMMModels[[3]]["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev




