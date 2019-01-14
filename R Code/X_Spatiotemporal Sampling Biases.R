
# Looking at spatiotemporal sampling biases ####

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
                             nchar(AssocsBase2$Year))%in%c("a","b"),"Year"],
         2,5)


AssocsBase2[substr(AssocsBase2$Year, 
                   nchar(AssocsBase2$Year),
                   nchar(AssocsBase2$Year))%in%c("a","b"),"Year"] <-
  substr(AssocsBase2[substr(AssocsBase2$Year, 
                            nchar(AssocsBase2$Year),
                            nchar(AssocsBase2$Year))%in%c("a","b"),"Year"],1,4)

AssocsBase2$Year <- as.numeric(AssocsBase2$Year)

AssocsDiscovery <- merge(AssocsBase2, 
                         SpatialHosts,
                         by.x = "Host", by.y = "Sp", all.x = T)

AssocsDiscovery <- merge(AssocsDiscovery, 
                         SpatialViruses,
                         by.x = "Virus", by.y = "Sp", all.x = T,
                         suffixes = c(".Host", ".Virus"))

ggplot(data.frame(x = as.numeric(names(table(AssocsDiscovery$Year))),
                  y = cumsum(table(AssocsDiscovery$Year))),
       aes(x, y)) +
  geom_point() + geom_smooth() +
  labs(x = "Year", y = "Discoveries")

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

ggplot(AssocsDiscovery, aes(LongMean.Host, LatMean.Host)) + 
  #geom_path(data = WorldMap, aes(long, lat, group = group)) +
  facet_wrap(~Year) + 
  theme(legend.position = "none") +
  geom_point(size = 5, aes(colour = as.factor(Human)), alpha = 0.1) +
  scale_colour_manual(values = colfunc(2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Year", 
       title = "Viral Discovery Host Location in Space and Time") +
  coord_fixed()

ggplot(AssocsDiscovery, aes(LongMean.Virus, LatMean.Virus)) + 
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
          "Human", "Domestic",
          "Year")

TestAssocs <- AssocsDiscovery[!NARows(AssocsDiscovery[,Vars]),]
TestAssocs <- TestAssocs %>% filter(Year>1955)

TestAssocs[,c("LongMean.Host","LatMean.Host")] <- TestAssocs[,c("LongMean.Host","LatMean.Host")]/50000

HostLocations = cbind(TestAssocs$LongMean.Host, TestAssocs$LatMean.Host)

WorldMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)

A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations) # Making A matrix
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

DiscoveryStack <- inla.stack(
  data = list(y = TestAssocs[,c("Human")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.index)) # Leave

DiscoveryIM[[1]] <- inla(f1, # Base model (no random effects)
                         family = c("binomial"),
                         data = inla.stack.data(DiscoveryStack),
                         control.compute = list(dic = TRUE),
                         control.predictor = list(A = inla.stack.A(DiscoveryStack))
)

DiscoveryIM[[2]] <- inla(f2, # f1 + SPDE random effect 
                         family = c("binomial"),
                         data = inla.stack.data(DiscoveryStack),
                         control.compute = list(dic = TRUE),
                         control.predictor = list(A = inla.stack.A(DiscoveryStack))
)

ggField(DiscoveryIM[[2]], WorldMesh) + 
  scale_fill_brewer(palette = AlberPalettes[3]) +
  geom_path(data = WorldMap/50000, inherit.aes = F, aes(long, lat, group = group)) +
  geom_point(data = TestAssocs, aes(LongMean.Host, LatMean.Host), inherit.aes = F) +
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
A3 <- inla.spde.make.A(WorldMesh, loc = HostLocations,
                       repl = as.numeric(TestAssocs$Group),
                       n.repl = NGroup) # Making A matrix

w.st <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.repl = NGroup)  

DiscoveryStack2 <- inla.stack(
  data = list(y = TestAssocs[,c("Human")]),  
  A = list(1, A3), # Vector of Multiplication factors              
  effects = list(
    Intercept = rep(1, N), # Leave
    #X = X, # Leave
    w = w.st)) # Leave

DiscoveryIM[[3]] <- inla(f3, # f1 + annual SPDE random effect w/o correlations 
                         family = c("binomial"),
                         data = inla.stack.data(DiscoveryStack2),
                         control.compute = list(dic = TRUE),
                         control.predictor = list(A = inla.stack.A(DiscoveryStack2))
)

DiscoveryIM[[4]] <- inla(f4, # f1 + annual SPDE random effect 
                         family = c("binomial"),
                         data = inla.stack.data(DiscoveryStack2),
                         control.compute = list(dic = TRUE),
                         control.predictor = list(A = inla.stack.A(DiscoveryStack2))
)


ggField(DiscoveryIM[[2]], WorldMesh, Groups = NGroups) + 
  scale_fill_brewer(palette = AlberPalettes[3]) +
  geom_path(data = WorldMap/50000, inherit.aes = F,
            aes(long, lat))

