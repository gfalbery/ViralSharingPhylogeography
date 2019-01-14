
# Preliminary spatial analysis ####

# VERY ROUGH, mostly copied across from deer analyses ####

library(INLA)

TestHosts <- SpatialHosts

TestHosts$hEigenvector <- kader:::cuberoot(TestHosts$Eigenvector)

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000

HCentSel <- INLAModelSel(SpatialHosts, "Eigenvector")

HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)

WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)

plot(WorldMesh)

A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix

spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1")), #paste(lapply(INLASelectList, function(b) b$Removed[[length(b$Removed)]])[[a]], collapse=" + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X<-as.data.frame(Xm[,]) # Model Matrix

#f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
#f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid')"))
#f3 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid') + f(w, model = spde)"))

f1 <- as.formula("y ~ -1 + Intercept" )
f2 <- as.formula(paste0("y ~ -1 + Intercept + ","f(Name, model = 'iid') + f(fYear, model= 'iid')"))
f3 <- as.formula(paste0("y ~ -1 + Intercept + ", "f(w, model = spde)"))

a = 2

EHAStack <- inla.stack(
  data = list(y = TestHosts[,c("Degree", "hEigenvector", "GroupSize")[a]]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    #Name=TestHosts$Name,
    #fYear=TestHosts$fYear,
    w = w.index)) # Leave

I1<-inla(f1, # Base model (no random effects)
         family=c("gaussian", "gaussian", "gaussian")[a],
         data = inla.stack.data(EHAStack),
         control.compute = list(dic = TRUE),
         control.predictor = list(A = inla.stack.A(EHAStack))
)

I3<-inla(f3, # f2 + SPDE random effect 
         family=c("gaussian", "gaussian", "gaussian")[a],
         data = inla.stack.data(EHAStack),
         control.compute = list(dic = TRUE),
         control.predictor = list(A = inla.stack.A(EHAStack))
)

# Plotting ####

ggField(I3, WorldMesh) + scale_fill_brewer(palette = AlberPalettes[1]) +
  geom_path(data = WorldMap/50000, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = TestHosts, aes(LongMean, LatMean), inherit.aes = F)

# Virus Spatial Analysis ####

TestViruses <- SpatialViruses

TestViruses$vEigenvector <- kader:::cuberoot(TestViruses$Eigenvector)

TestViruses[,c("LongMean","LatMean")] <- TestViruses[,c("LongMean","LatMean")]/50000

VirusLocations = cbind(TestViruses$LongMean, TestViruses$LatMean)

A3 <- inla.spde.make.A(WorldMesh, loc = VirusLocations) # Making A matrix

spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1")), #paste(lapply(INLASelectList, function(b) b$Removed[[length(b$Removed)]])[[a]], collapse=" + "))), 
                   data = TestViruses)
N <- nrow(TestViruses)
X<-as.data.frame(Xm[,]) # Model Matrix

#f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
#f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid')"))
#f3 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid') + f(w, model = spde)"))

f1 <- as.formula("y ~ -1 + Intercept" )
f2 <- as.formula(paste0("y ~ -1 + Intercept + ","f(Name, model = 'iid') + f(fYear, model= 'iid')"))
f3 <- as.formula(paste0("y ~ -1 + Intercept + ", "f(w, model = spde)"))

a = 2

ViralStack <- inla.stack(
  data = list(y = TestViruses[,c("Degree", "vEigenvector", "GroupSize")[a]]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    #Name=TestViruses$Name,
    #fYear=TestViruses$fYear,
    w = w.index)) # Leave

I1<-inla(f1, # Base model (no random effects)
         family=c("gaussian", "gaussian", "gaussian")[a],
         data = inla.stack.data(ViralStack),
         control.compute = list(dic = TRUE),
         control.predictor = list(A = inla.stack.A(ViralStack))
)

I3<-inla(f3, # f2 + SPDE random effect 
         family=c("gaussian", "gaussian", "gaussian")[a],
         data = inla.stack.data(ViralStack),
         control.compute = list(dic = TRUE),
         control.predictor = list(A = inla.stack.A(ViralStack))
)

# Plotting ####

ggField(I3, WorldMesh) + scale_fill_brewer(palette = AlberPalettes[1]) +
  geom_point(data = TestViruses, aes(LongMean, LatMean), inherit.aes = F)

# Either centroids don't do a great job of representing spatial distribution 
# of viruses, or the network really isn't that spatially autocorrelated.




