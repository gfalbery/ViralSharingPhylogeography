
# Spatial autocorrelation in space: Citations are highly autocorrelated ####

Resp = "hAllZACites"

CiteCovar1 = c(
  "hDom",
  "hOrder"
)

TestHosts <- Hosts %>% dplyr::select(Resp, CiteCovar1, "LongMean", "LatMean") %>%
  
  mutate(hAllZACites = log(hAllZACites + 1)
         #hDiseaseZACites = log(hDiseaseZACites + 1),
  ) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,CiteCovar1, "LongMean", "LatMean")])))

CiteModelList <- list()

# SPDE Model

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(CiteCovar1, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("hDomdomestic"))]) # Model Matrix

f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "),
                        " + f(w, model = spde)"))

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    w = w.index)) # Leave

CiteModelList[[1]] <- inla(f1, # Base model (no random effects)
                           family = c("nbinomial"),
                           data = inla.stack.data(CentStack),
                           control.compute = list(dic = TRUE),
                           control.predictor = list(A = inla.stack.A(CentStack))
)


CiteModelList[[2]] <- inla(f2, # Base model (no random effects)
                           family = c("nbinomial"),
                           data = inla.stack.data(CentStack),
                           control.compute = list(dic = TRUE),
                           control.predictor = list(A = inla.stack.A(CentStack))
)

sapply(CiteModelList, function(a) a$dic$dic) %>% diff

ggField(CiteModelList[[2]], WorldMesh) + 
  geom_path(data = WorldMap[,c("long", "lat", "group")]/50000, inherit.aes = F, aes(long, lat, group = group)) + 
  geom_point(data = TestHosts, inherit.aes = F, aes(x = LongMean, y = LatMean)) +
  scale_fill_brewer(palette = AlberPalettes[2])

# Spatial autocorrelation in space ####

Resp = "hDiseaseZACites"

CiteCovar1 = c(
  "hDom",
  "hOrder"
)

TestHosts <- Hosts %>% dplyr::select(Resp, CiteCovar1, "LongMean", "LatMean") %>%
  
  mutate(#hAllZACites = log(hAllZACites + 1),
         hDiseaseZACites = log(hDiseaseZACites + 1)
  ) %>%
  
  slice(which(!NARows(Hosts[,c(Resp,CiteCovar1, "LongMean", "LatMean")])))

# SPDE Model

TestHosts[,c("LongMean","LatMean")] <- TestHosts[,c("LongMean","LatMean")]/50000
HostLocations = cbind(TestHosts$LongMean, TestHosts$LatMean)
WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(CiteCovar1, collapse = " + "))), 
                   data = TestHosts)
N <- nrow(TestHosts)
X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("hDomdomestic"))]) # Model Matrix

f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "),
                        " + f(w, model = spde)"))

CentStack <- inla.stack(
  data = list(y = TestHosts[,Resp]),  
  A = list(1, 1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    w = w.index)) # Leave

CiteModelList[[3]] <- inla(f1, # Base model (no random effects)
                           family = c("nbinomial"),
                           data = inla.stack.data(CentStack),
                           control.compute = list(dic = TRUE),
                           control.predictor = list(A = inla.stack.A(CentStack))
)


CiteModelList[[4]] <- inla(f2, # Base model (no random effects)
                           family = c("nbinomial"),
                           data = inla.stack.data(CentStack),
                           control.compute = list(dic = TRUE),
                           control.predictor = list(A = inla.stack.A(CentStack))
)

sapply(CiteModelList, function(a) a$dic$dic) %>% diff

ggField(CiteModelList[[4]], WorldMesh) + 
  scale_fill_brewer(palette = AlberPalettes[2]) +
  geom_path(data = WorldMap[,c("long", "lat", "group")]/50000, inherit.aes = F, aes(long, lat, group = group)) + 
  geom_point(data = TestHosts, inherit.aes = F, aes(x = LongMean, y = LatMean))

# Do they correlate?

qplot(CiteModelList[[2]]$summary.random$w$mean,CiteModelList[[4]]$summary.random$w$mean)

# Hells ye

SaveCiteModelList <- CiteModelList