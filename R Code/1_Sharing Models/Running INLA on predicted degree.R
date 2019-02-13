
# Running INLA on predicted degree

library(INLA); library(tidyverse); library(ggregplot)

TestHosts <- na.omit(Panth1[,c("AllPredDegree","Sp","LongMean","LatMean")]) %>%
  mutate(AllPredDegree = scale(AllPredDegree))

N = nrow(TestHosts)

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(50, 50), cutoff = 20)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Establishing model formulae ####
f1 = as.formula(paste("y ~ -1 + Intercept"))

f2 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        #paste(names(X), collapse = " + "),
                        "f(w, model = spde)"))

AllPredStack <- inla.stack(
  data = list(y = TestHosts$AllPredDegree),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

IM1 <-
  inla(
    f1, 
    family = "gaussian",
    data = inla.stack.data(AllPredStack),
    control.compute = list(dic = TRUE),
    control.predictor = list(A = inla.stack.A(AllPredStack))
  )

IM2 <-
  inla(
    f2, 
    family = "gaussian",
    data = inla.stack.data(AllPredStack),
    control.compute = list(dic = TRUE),
    control.predictor = list(A = inla.stack.A(AllPredStack))
  )

INLADICFig(list(IM1, IM2))

ggField(IM2, WorldMesh)

ggField(IM2, WorldMesh) + scale_fill_brewer(palette = AlberPalettes[2]) +
  geom_path(data = WorldMap, inherit.aes = F, aes(long/50000, lat/50000, group = group)) +
  geom_point(data = TestHosts, inherit.aes = F, aes(LongMean, LatMean), alpha = 0.2)

