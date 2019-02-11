
# Running INLA on predicted degree

library(INLA)
TestHosts <- Panth1[,c("AllPredDegree","MSW05_Binomial")]
N = nrow(TestHosts)

TestHosts$Sp <- TestHosts$MSW05_Binomial

Long <- with(FullPolygons, tapply(long, Host, mean))
Lat <- with(FullPolygons, tapply(lat, Host, mean))

TestHosts[,c("LongMean","LatMean")] <- cbind(Long[TestHosts$Sp],Lat[TestHosts$Sp])/50000
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

