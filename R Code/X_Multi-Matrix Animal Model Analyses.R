
# Fitting Animal Models using matrices ####

load("Model Files/HostMMModels.Rdata")

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

INLARep(HostMMModels[[2]], "IndexSpace")
INLARep(HostMMModels[[3]], "IndexPhylo")

# From joint models ####

# Extracting variance readings ####

MySqrt <- function(x) {1 / sqrt(x) }
tau <- HostMMModels[[4]]$marginals.hyperpar$`Precision for the Gaussian observations`
sigma <- inla.emarginal(MySqrt, tau)
sigma

sigma2 <- inla.emarginal(MySqrt, HostMMModels[[4]]$marginals.hyperpar[[paste0("Precision for ", "IndexSpace")]])
sigma3 <- inla.emarginal(MySqrt, HostMMModels[[4]]$marginals.hyperpar[[paste0("Precision for ", "IndexPhylo")]])

allvar <- sigma^2 + sigma2^2 + sigma3^2

Sigma1 <- sigma^2/allvar
Sigma2 <- sigma2^2/allvar
Sigma3 <- sigma3^2/allvar

list(IndexSpace = sigma2^2/(sigma^2 + sigma2^2 + sigma3^2),
     IndexPhylo = sigma3^2/(sigma^2 + sigma2^2 + sigma3^2))

df = data.frame(Model = c("Space", "Space", "Phylo", "Phylo", "Both", "Both", "Both"),
                Var = c("Resid", "Space", "Resid", "Phylo", "Resid", "Space", "Phylo"),
                Value = c(1-INLARep(HostMMModels[[2]], "IndexSpace")$Estimate, 
                          INLARep(HostMMModels[[2]], "IndexSpace")$Estimate,
                          1-INLARep(HostMMModels[[3]], "IndexPhylo")$Estimate,
                          INLARep(HostMMModels[[3]], "IndexPhylo")$Estimate,
                          Sigma1, Sigma2, Sigma3))

df$Var = factor(df$Var, levels = c("Resid", "Space", "Phylo"))
df$Model = factor(df$Model, levels = c("Space", "Phylo", "Both"))

ggplot(df, aes(Model, Value, fill = Var)) + geom_col(position = "stack") + 
  scale_fill_manual(values = c("grey", "#2C7FB8","#DE2D26")) + labs(y = "Variance accounted for", fill = "Component") +
  ggsave("Figures/Variance Components.jpeg", units = "mm", height = 100, width = 100, dpi = 300)

# Phylogenetic relatedness feeds more into centrality in the network than spatial sharing does???

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




