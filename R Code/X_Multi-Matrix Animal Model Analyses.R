
# Fitting Animal Models using matrices ####

# Making Host Data ####

SpaceContacts <- as(solve(RangeAdj1[FHN, FHN]),"dgCMatrix")
GRMatrix <- as(solve(1-CytBMatrix[FHN, FHN]),"dgCMatrix")

# NB phylo matrix has been inverted ####

FHosts$IndexSpace = sapply(FHosts$Sp, function(a) which(FHN==a))
FHosts$IndexPhylo = sapply(FHosts$Sp, function(a) which(FHN==a))

TestHosts <- FHosts

TestHosts$hEigenvector <- kader:::cuberoot(TestHosts$Eigenvector)

f1 = hEigenvector ~ 1

M1 <- inla(f1, 
           data = TestHosts, 
           family = "gaussian",
           control.compute = list(dic = T))


f2 =  hEigenvector ~ f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))

M2 <- inla(f2, 
           data = TestHosts, 
           family = "gaussian",
           control.compute = list(dic = T))


f3 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5))

M3 <- inla(f3, 
           data = TestHosts, 
           family = "gaussian",
           control.compute = list(dic = T))

f4 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

M4 <- inla(f4, 
           data = TestHosts, 
           family = "gaussian",
           control.compute = list(dic = T))

Models <- list(M1, M2, M3)#, M4)

sapply(Models, function(a) a$dic$dic)

# Fitting space sharing and phylogenies both improve centrality fit 

# Extracting variance readings ####

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  M2$marginals.hyperpar$"Precision for IndexSpace")

e.IndexA = inla.expectation(function(x) x, sigma.IndexSpace)
ci.IndexA = inla.qmarginal(c(0.025, 0.975), sigma.IndexSpace)

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  M2["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace = 1/inla.qmarginal(c(0.025, 0.975), 
                             M2["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo = inla.emarginal(function(x) 1/x,
                                  M3["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo = 1/inla.qmarginal(c(0.025, 0.975), 
                             M3["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev


# Making Virus dataset ####

VirusSpaceContacts <- as(solve(VirusRangeAdj1[FVN, FVN]+matrix(runif(length(FVN)^2, min=0,max=0.01), ncol = length(FVN))),"dgCMatrix")
#GRMatrix <- as(solve(1-CytBMatrix[FVN, FVN]),"dgCMatrix")

# NB phylo matrix has been inverted ####

FViruses$IndexSpace = sapply(FViruses$Sp, function(a) which(FVN==a))

TestViruses <- FViruses

TestViruses$vEigenvector <- kader:::cuberoot(TestViruses$Eigenvector)

f1 = vEigenvector ~ 1

M1 <- inla(f1, 
           data = TestViruses, 
           family = "gaussian",
           control.compute = list(dic = T))

f2 =  hEigenvector ~ f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))

M2 <- inla(f2, 
           data = TestViruses, 
           family = "gaussian",
           control.compute = list(dic = T))


f3 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5))

M3 <- inla(f3, 
           data = TestViruses, 
           family = "gaussian",
           control.compute = list(dic = T))

f4 = hEigenvector ~ f(IndexPhylo, model = "generic0", Cmatrix = GRMatrix,
                      constr = TRUE,param = c(0.5, 0.5)) + 
  f(IndexSpace, model="generic0", Cmatrix = VirusSpaceContacts,
    constr=TRUE,param = c(0.5, 0.5))

M4 <- inla(f4, 
           data = TestViruses, 
           family = "gaussian",
           control.compute = list(dic = T))

Models <- list(M1, M2)#, M3)#, M4)

sapply(Models, function(a) a$dic$dic)

# Fitting space sharing and phylogenies both improve centrality fit 

# Extracting variance readings ####

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  M2$marginals.hyperpar$"Precision for IndexSpace")

e.IndexA = inla.expectation(function(x) x, sigma.IndexSpace)
ci.IndexA = inla.qmarginal(c(0.025, 0.975), sigma.IndexSpace)

sigma.IndexSpace = inla.emarginal(function(x) 1/x,
                                  M2["marginals.hyperpar"][[1]][["Precision for IndexSpace"]])

ci.IndexSpace = 1/inla.qmarginal(c(0.025, 0.975), 
                                 M2["marginals.hyperpar"][[1]][["Precision for IndexSpace"]]) %>%
  rev

sigma.IndexPhylo = inla.emarginal(function(x) 1/x,
                                  M3["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]])

ci.IndexPhylo = 1/inla.qmarginal(c(0.025, 0.975), 
                                 M3["marginals.hyperpar"][[1]][["Precision for IndexPhylo"]]) %>%
  rev




