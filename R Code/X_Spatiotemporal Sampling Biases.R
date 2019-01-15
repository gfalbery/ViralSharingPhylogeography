
# Looking at spatiotemporal sampling biases ####

load("Model Files/Discovery Models.Rdata")

NARows <-function(df, vars){
  apply(as.data.frame(df[,vars]), 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

# Data frame replication ####

AssocsBase2 <- AssocsBase
AssocsBase2 <- droplevels(AssocsBase[!AssocsBase$Host == "Homo_sapiens"&
                                       !AssocsBase$Virus == "Rabies_virus",])

AssocsBase2$Year <- as.character(substr(AssocsBase2$Reference, 
                                        nchar(AssocsBase2$Reference)-4, 
                                        nchar(AssocsBase2$Reference)))

AssocsBase2[!substr(AssocsBase2$Year, 
                    nchar(AssocsBase2$Year),
                    nchar(AssocsBase2$Year))%in%c("a","b"),"Year"] <-
  substr(AssocsBase2[!substr(AssocsBase2$Year, 
                             nchar(AssocsBase2$Year),
                             nchar(AssocsBase2$Year))%in%c("a","b"),"Year"],2,5)

AssocsBase2[substr(AssocsBase2$Year, 
                   nchar(AssocsBase2$Year),
                   nchar(AssocsBase2$Year))%in%c("a","b"),"Year"] <-
  substr(AssocsBase2[substr(AssocsBase2$Year, 
                            nchar(AssocsBase2$Year),
                            nchar(AssocsBase2$Year))%in%c("a","b"),"Year"],1,4)

AssocsBase2$Year <- as.numeric(AssocsBase2$Year)

AssocsDiscovery <- merge(AssocsBase2, 
                         Hosts,
                         by.x = "Host", by.y = "Sp", all.x = T)

AssocsDiscovery <- merge(AssocsDiscovery, 
                         Viruses,
                         by.x = "Virus", by.y = "Sp", all.x = T,
                         suffixes = c(".Host", ".Virus"))

AssocsDiscovery <- AssocsDiscovery %>% dplyr::rename(LongMean.Host = LongMean, LatMean.Host = LatMean)

ggplot(data.frame(x = as.numeric(names(table(AssocsDiscovery$Year))),
                  y = cumsum(table(AssocsDiscovery$Year))),
       aes(x, y)) +
  geom_point() + geom_smooth(colour = AlberColours[3]) +
  labs(x = "Year", y = "HP3 References", title = "Viral Reports Over Time") +
  ggsave("Figures/Reference discovery over time.jpeg", 
         units = "mm", width= 100, height = 100, dpi = 300)  

colfunc <- colorRampPalette(c("blue", "red"))

ggplot(AssocsDiscovery, aes(LongMean.Host, LatMean.Host)) + 
  geom_path(data = WorldMap, aes(long, lat, group = group)) +
  #facet_wrap(~Year) + 
  #theme(legend.position = "none") +
  geom_point(size = 5, aes(colour = as.factor(Year)), alpha = 0.1) +
  scale_colour_manual(values = colfunc(nunique(AssocsDiscovery$Year))) +
  labs(x = "Longitude", y = "Latitude", colour = "Year", 
       title = "Viral Discovery Host Location in Space and Time") +
  coord_fixed() +
  ggsave("Figures/Hosts in Space and Time.jpeg", 
         units = "mm", width= 200, height = 150, dpi = 300)

ggplot(AssocsDiscovery, aes(LongMean.Total, LatMean.Total)) + 
  geom_path(data = WorldMap, aes(long, lat, group = group)) +
  #facet_wrap(~Year) + 
  #theme(legend.position = "none") +
  geom_point(size = 5, aes(colour = as.factor(Year)), alpha = 0.1) +
  scale_colour_manual(values = colfunc(nunique(AssocsDiscovery$Year))) +
  labs(x = "Longitude", y = "Latitude", colour = "Year", 
       title = "Viral Discovery Location in Space and Time") +
  coord_fixed() +
  ggsave("Figures/Viral Discovery Location in Space and Time.jpeg", 
         units = "mm", width= 200, height = 150, dpi = 300)

# INLA Models: zoonosis probability ####

library(INLA)

Vars <- c("LongMean.Host", "LatMean.Host",
          #"LongMean.Centroid", "LatMean.Centroid",
          "Human", "Domestic.Virus",
          "Year")

TestAssocs <- AssocsDiscovery[!NARows(AssocsDiscovery[,Vars]),]
TestAssocs <- TestAssocs %>% filter(Year>1955)

TestAssocs[,c("LongMean.Host","LatMean.Host")] <- TestAssocs[,c("LongMean.Host","LatMean.Host")]/50000
# TestAssocs[,c("LongMean.Centroid","LatMean.Centroid")] <- TestAssocs[,c("LongMean.Centroid","LatMean.Centroid")]/50000

TestAssocs$Group <- cut(TestAssocs$Year, breaks = rep(0:1000)*5, labels = 1:1000) %>% 
  droplevels %>% factor %>% as.numeric
TestAssocs$Group[TestAssocs$Group==13] <- 12
NGroup <- nunique(TestAssocs$Group)

Locations = cbind(TestAssocs$LongMean.Host, TestAssocs$LatMean.Host)
# HostLocations = cbind(TestAssocs$LongMean.Centroid, TestAssocs$LatMean.Centroid)

WorldMesh <- inla.mesh.2d(loc = Locations, max.edge = c(10, 25), cutoff = 10)

A3 <- inla.spde.make.A(WorldMesh, loc = Locations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1")), #paste(lapply(INLASelectList, function(b) b$Removed[[length(b$Removed)]])[[a]], collapse=" + "))), 
                   data = TestAssocs)
N <- nrow(TestAssocs)
X<-as.data.frame(Xm[,]) # Model Matrix

#f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
#f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid')"))
#f3 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(Name, model = 'iid') + f(fYear, model= 'iid') + f(w, model = spde)"))

f1 <- as.formula("y ~ -1 + Intercept" )
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", "f(w, model = spde)"))

ZooDiscoveryStack <- inla.stack(
  data = list(y = TestAssocs[,c("Human")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

ZooDiscoveryList <- list()

ZooDiscoveryList[[1]] <- inla(f1, # Base model (no random effects)
                              family = c("binomial"),
                              data = inla.stack.data(ZooDiscoveryStack),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(ZooDiscoveryStack))
)

ZooDiscoveryList[[2]] <- inla(f2, # f1 + SPDE random effect 
                              family = c("binomial"),
                              data = inla.stack.data(ZooDiscoveryStack),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(ZooDiscoveryStack))
)

ggField(ZooDiscoveryList[[2]], WorldMesh) + 
  scale_fill_brewer(palette = AlberPalettes[3]) +
  geom_path(data = WorldMap/50000,  inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = dplyr::select(TestAssocs, -Group), aes(LongMean.Host, LatMean.Host), inherit.aes = F) +
  labs(x  = "Longitude", y = "Latitude", fill = "Zoonosis", 
       title = "Spatial Autocorrelation in Human Infection Probability") +
  ggsave("Figures/Human Virus Discovery INLA Effect.jpeg", 
         units = "mm", height = 150, width = 200, dpi = 300)

# Trying spatiotemporal model ####

f3 <- y ~ -1 + Intercept + 
  f(w, model = spde, replicate = w.repl)

f4 <- y ~ -1 + Intercept + 
  f(w, model = spde, 
    group = w.group,
    control.group = list(model="iid"))

TestAssocs$Group <- cut(TestAssocs$Year, breaks = rep(0:1000)*5, labels = 1:1000) %>% 
  droplevels %>% factor %>% as.numeric
TestAssocs$Group[TestAssocs$Group==13] <- 12
NGroup <- nunique(TestAssocs$Group)

spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
A3 <- inla.spde.make.A(WorldMesh, loc = Locations,
                       repl = as.numeric(TestAssocs$Group),
                       n.repl = NGroup) # Making A matrix

w.st <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.repl = NGroup)  

ZooDiscoveryStack2 <- inla.stack(
  data = list(y = TestAssocs[,c("Human")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.st)) # Leave

ZooDiscoveryList[[3]] <- inla(f3, # f1 + annual SPDE random effect w/o correlations 
                              family = c("binomial"),
                              data = inla.stack.data(ZooDiscoveryStack2),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(ZooDiscoveryStack2))
)

ZooDiscoveryList[[4]] <- inla(f4, # f1 + annual SPDE random effect 
                              family = c("binomial"),
                              data = inla.stack.data(ZooDiscoveryStack2),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(ZooDiscoveryStack2))
)

Labels <- as.character(paste0(1:NGroup*5 +1950, "-", (1:NGroup+1)*5 +1950))
names(Labels) <- 1:NGroup

ggField(ZooDiscoveryList[[3]], WorldMesh, Groups = NGroup) + 
  facet_wrap(~Group, labeller = labeller(Group = Labels), nrow = 4) +
  scale_fill_brewer(palette = AlberPalettes[3]) +
  geom_path(data = WorldMap/50000,  colour = "black",inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = TestAssocs, aes(LongMean.Host, LatMean.Host), inherit.aes = F) +
  labs(x  = "Longitude", y = "Latitude", fill = "Zoonosis", 
       title = "Spatiotemporal Autocorrelation in Human Infection Probability") +
  ggsave("Figures/Spatiotemporal Zoonosis Discovery INLA Effect.jpeg", 
         units = "mm", height = 300, width = 350, dpi = 300)

sapply(ZooDiscoveryList, function(a) a$dic$dic)

# INLA Models: domestic infection probability ####

A3 <- inla.spde.make.A(WorldMesh, loc = Locations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

DomDiscoveryStack <- inla.stack(
  data = list(y = TestAssocs[,c("Domestic.Virus")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

DomDiscoveryList <- list()

DomDiscoveryList[[1]] <- inla(f1, # Base model (no random effects)
                              family = c("binomial"),
                              data = inla.stack.data(DomDiscoveryStack),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(DomDiscoveryStack))
)

DomDiscoveryList[[2]] <- inla(f2, # f1 + SPDE random effect 
                              family = c("binomial"),
                              data = inla.stack.data(DomDiscoveryStack),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(DomDiscoveryStack))
)

sapply(DomDiscoveryList, function(a) a$dic$dic)

ggField(DomDiscoveryList[[2]], WorldMesh) + 
  scale_fill_brewer(palette = AlberPalettes[3]) +
  geom_path(data = WorldMap/50000,  inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = dplyr::select(TestAssocs, -Group), aes(LongMean.Host, LatMean.Host), inherit.aes = F) +
  labs(x  = "Longitude", y = "Latitude", fill = "Domestic", 
       title = "Spatial Autocorrelation in Domestic Virus Discovery") +
  ggsave("Figures/Domestic Virus Discovery INLA Effect.jpeg", 
         units = "mm", height = 150, width = 200, dpi = 300)

# Trying spatiotemporal model ####

f3 <- y ~ -1 + Intercept + 
  f(w, model = spde, replicate = w.repl)

f4 <- y ~ -1 + Intercept + 
  f(w, model = spde, 
    group = w.group,
    control.group = list(model="iid"))

spde = inla.spde2.pcmatern(mesh = WorldMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
A3 <- inla.spde.make.A(WorldMesh, loc = Locations,
                       repl = as.numeric(TestAssocs$Group),
                       n.repl = NGroup) # Making A matrix

w.st <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.repl = NGroup)  

DomDiscoveryStack2 <- inla.stack(
  data = list(y = TestAssocs[,c("Domestic.Virus")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.st)) # Leave

DomDiscoveryList[[3]] <- inla(f3, # f1 + annual SPDE random effect w/o correlations 
                              family = c("binomial"),
                              data = inla.stack.data(DomDiscoveryStack2),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(DomDiscoveryStack2))
)

DomDiscoveryList[[4]] <- inla(f4, # f1 + annual SPDE random effect 
                              family = c("binomial"),
                              data = inla.stack.data(DomDiscoveryStack2),
                              control.compute = list(dic = TRUE),
                              control.predictor = list(A = inla.stack.A(DomDiscoveryStack2))
)

ggField(DomDiscoveryList[[3]], WorldMesh, Groups = NGroup) + 
  facet_wrap(~Group, labeller = labeller(Group = Labels), nrow = 4) +
  scale_fill_brewer(palette = AlberPalettes[3]) +
  geom_path(data = WorldMap/50000,  inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = TestAssocs, aes(LongMean.Host, LatMean.Host), inherit.aes = F) +
  labs(x  = "Longitude", y = "Latitude", fill = "Domestic", 
       title = "Spatial Autocorrelation in Domestic Infection Probability") +
  ggsave("Figures/Spatiotemporal Domestic Virus Discovery INLA Effect.jpeg", 
         units = "mm", height = 200, width = 300, dpi = 300)

sapply(DomDiscoveryList, function(a) a$dic$dic)

DiscoveryModelList <- list(ZooDiscoveryList, DomDiscoveryList)

save(ZooDiscoveryList, file = "Model Files/Human Discovery Models.Rdata")
save(DomDiscoveryList, file = "Model Files/Domestic Discovery Models.Rdata")
