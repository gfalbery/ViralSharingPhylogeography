
# Stuff to show Noam

# summary()


# Terrible chain convergence ####

i = 1

mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$Sol[,1:14])
#mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$VCV)
mc <- do.call(mcmc.list, mc)
#par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
#gelman.plot(mc, auto.layout=F)
#gelman.diag(mc)
par(mfrow=c(7,4), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

# Only slightly better w/o G ####

i = 2

mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$Sol[,1:14])
#mc <- ZI_runs[1:10 + 10*(i-1)] %>% lapply(function(a) a$VCV)
mc <- do.call(mcmc.list, mc)
#par(mfrow=c(7,2), mar=c(2,2,1,2), ask = F)
#gelman.plot(mc, auto.layout=F)
#gelman.diag(mc)
par(mfrow=c(7,4), mar=c(2, 1, 1, 1))
plot(mc, ask=F, auto.layout=F)

# Effects largely in agreement within model ####

Efxplot(ZI_runs[1:10])

# But massive differences between models ####

Efxplot(ModelList)

ggplot(ZISols, aes(x = Variable, y = Estimate, colour = as.factor(Model))) + 
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(position = position_dodge(w = 0.5), aes(ymin = Lower, 
                                                        ymax = Upper), size = 0.3, width = 0.2) + 
  geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL, colour = "Model") + coord_flip() + AlberTheme +
  facet_wrap(~Component, scales = "free_x")

# Fit of model 1 predictions ####

ggplot(FinalHostMatrix, aes(Virus, PredVirus1)) + geom_point() + geom_smooth() + coord_fixed()

ggplot(Hosts, aes(Degree, PredDegree1)) + geom_point() + geom_smooth() + coord_fixed()

# Fit of model 1 predictions sans random ####

ggplot(FinalHostMatrix, aes(Virus, PredVirus1b)) + geom_point() + geom_smooth()+ coord_fixed()

ggplot(Hosts, aes(Degree, PredDegree1b)) + geom_point() + geom_smooth() + coord_fixed()

# Fit of model 2 predictions #####

ggplot(FinalHostMatrix, aes(Virus, PredVirus2)) + geom_point() + geom_smooth() + coord_fixed()

ggplot(Hosts, aes(Degree, PredDegree2)) + geom_point() + geom_smooth() + coord_fixed()

# Fit of model 1 predictions with increased citations ####

ggplot(FinalHostMatrix, aes(Virus, PredVirus3)) + geom_point() + geom_smooth() + coord_fixed()

ggplot(Hosts, aes(Degree, PredDegree3)) + geom_point() + geom_smooth() + coord_fixed()

# Pairs of virals ####

GGally::ggpairs(FinalHostMatrix[,c("Virus","PredVirus1","PredVirus1b","PredVirus2","PredVirus3")], 
                lower = list(continuous = "smooth"))

# Pairs of degrees ####

GGally::ggpairs(Hosts[,c("Degree","PredDegree1","PredDegree1b","PredDegree2","PredDegree3")], 
                lower = list(continuous = "smooth"))



