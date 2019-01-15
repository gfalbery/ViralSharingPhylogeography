
# Preliminary spatial analysis ####

# VERY ROUGH, mostly copied across from deer analyses ####

library(INLA)

HostCentralityList <- list()

TestHosts <- SpatialHosts

TestHosts$hEigenvector <- kader:::cuberoot(TestHosts$Eigenvector)

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000

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

f1 <- as.formula("y ~ -1 + Intercept" )
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", "f(w, model = spde)"))

a = 2

EHAStack <- inla.stack(
  data = list(y = TestHosts[,c("hEigenvector")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

HostCentralityList[[1]] <- inla(f1, # Base model (no random effects)
                                family = c("gaussian"),
                                data = inla.stack.data(EHAStack),
                                control.compute = list(dic = TRUE),
                                control.predictor = list(A = inla.stack.A(EHAStack))
)

HostCentralityList[[2]] <- inla(f2, # f2 + SPDE random effect 
                                family = "gaussian",
                                data = inla.stack.data(EHAStack),
                                control.compute = list(dic = TRUE),
                                control.predictor = list(A = inla.stack.A(EHAStack))
)

# Plotting ####

ggField(HostCentralityList[[2]], WorldMesh) + scale_fill_brewer(palette = AlberPalettes[1]) +
  geom_path(data = WorldMap/50000, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = TestHosts, aes(LongMean, LatMean), inherit.aes = F) +
  labs(x  = "Longitude", y = "Latitude", fill = "Centrality", 
       title = "Spatial Autocorrelation in Host Centrality") +
  ggsave("Figures/Host Centrality INLA Effect.jpeg", 
         units = "mm", height = 150, width = 200, dpi = 300)

# Virus Spatial Analysis ####

VirusCentralityList <- list()

TestViruses <- SpatialViruses

TestViruses$vEigenvector <- kader:::cuberoot(TestViruses$Eigenvector)

TestViruses[,c("LongMean","LatMean")] <- TestViruses[,c("LongMean.Total","LatMean.Total")]/50000

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
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", "f(w, model = spde)"))

a = 2

ViralStack <- inla.stack(
  data = list(y = TestViruses[,c("vEigenvector")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

VirusCentralityList[[1]] <- inla(f1, # Base model (no random effects)
                                 family = "gaussian",
                                 data = inla.stack.data(ViralStack),
                                 control.compute = list(dic = TRUE),
                                 control.predictor = list(A = inla.stack.A(ViralStack))
)

VirusCentralityList[[2]] <- inla(f2, # f2 + SPDE random effect 
                            family = "gaussian",
                            data = inla.stack.data(ViralStack),
                            control.compute = list(dic = TRUE),
                            control.predictor = list(A = inla.stack.A(ViralStack))
)

sapply(VirusCentralityList, function(a) a$dic$dic) # Slight Improvement

# Plotting ####

ggField(VirusCentralityList[[2]], WorldMesh) + scale_fill_brewer(palette = AlberPalettes[2]) +
  geom_path(data = WorldMap/50000, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = TestViruses, aes(LongMean, LatMean), inherit.aes = F) +
  labs(x  = "Longitude", y = "Latitude", fill = "Centrality", 
       title = "Weak Spatial Autocorrelation in Virus Centrality") +
  ggsave("Figures/Virus Centrality INLA Effect.jpeg", 
         units = "mm", height = 150, width = 200, dpi = 300)

# Either centroids don't do a great job of representing spatial distribution 
# of viruses, or the network really isn't that spatially autocorrelated.




