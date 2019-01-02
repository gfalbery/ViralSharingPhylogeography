# Elaborating on data frames ####

# Which viruses infect domestic hosts? ####

dim(Viruses)

library(magrittr); library(dplyr); library(RColorBrewer)

Viruses <- Viruses %>%
  mutate(
    Domestic = case_when(
      Sp %in% DomesticViruses ~ 1,
      TRUE ~ 0),
    
    PropDomestic = c(table(AssocsTraits[AssocsTraits$Domestic == 1, "Virus"])/
                       table(AssocsTraits$Virus)),
    
    Records = c(table(AssocsTraits$Virus))
  )

ggMMplot(Viruses, "vFamily", "Domestic")
ggMMplot(Viruses, "Human", "Domestic")

ggplot(Viruses, aes(PropDomestic, Human)) + 
  geom_point(colour = AlberColours[5], alpha = 0.3) + 
  geom_smooth(colour = "black", fill= NA) + 
  stat_smooth(colour = "black", geom = "ribbon", fill = NA, lty = 2) +
  #geom_smooth(method = lm, formula = y ~ poly(x,2), colour = "black", fill = NA) +
  #stat_smooth(method = lm, formula = y ~ poly(x,2), colour = "black", geom = "ribbon", fill = NA, lty = 2) +
  coord_fixed() + 
  labs(x = "Proportion Domestic Hosts", y = "Zoonotic", title = "Are domestic viruses more likely to be zoonotic?") +
  theme(plot.title = element_text(hjust = 0.5))

# Seems like a negative trend - weird. The hump in the middle may be the result of 
# high host specificity of those with low proportion domestic hosts?

# Testing dom-zoo associations ####
# Does ```records``` correlate with host specificity? 

resp1 <- "IsZoonotic"
resp2 <- "PropDomestic"

Viruscovar <- c("Records","vCytoReplicTF","vDNAoRNA","vEnvelope","vGenomeAveLengthLn","vPubMedCitesLn","vSegmentedTF","vSSoDS","vVectorYNna", "PropDomestic")
Hostcovar <- c("UrbRurPopRatioLn","AreaHost","HabInhabitedChgLn","hAllZACitesLn","hMarOTerr","hMassGramsPVR","HumPopDensLnChg","TotHumPopLn")
Domcovar <- c("domestic_category","DOMYearBP")

DomZooFormula1 <- as.formula(paste0(resp1, " ~ ", paste(c(Viruscovar,resp2), collapse = " + ")))
DomZooFormula2 <- as.formula(paste0("cbind(",resp1,",",resp2,")", " ~ ", paste(c(covar), collapse = " + ")))

DomZooPrior <- 
  list(R = list(V = diag(1), nu = 0.002, fix = 1))

DomZooMultivPrior <- 
  list(R = list(V = diag(2), nu = 0.002, fix = 2))

DomZooPrior2 <- 
  list(R = list(V = diag(1), nu = 0.002, fix = 1),
       G = list(G1 = list(V = diag(1), nu = 0.002)))

DomZooMultivPrior2 <- 
  list(R = list(V = diag(2), nu = 0.002, fix = 2),
       G = list(G1 = list(V = diag(2), nu = 0.002)))

mf <- 2 # Multiplication factor for MCMC iterations

DomZooMCMC1 <- MCMCglmm(
  DomZooFormula1,
  data = Viruses,
  prior = DomZooPrior,
  family = "categorical",
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
) 

DomZooMultiv1 <- MCMCglmm(
  cbind(PropDomestic, IsZoonotic) ~ trait - 1 + trait:(Records),
  prior = DomZooMultivPrior,
  rcov =~ us(trait):units,
  data = Viruses,
  family = c("gaussian","categorical"),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
) 

DomZooMCMC2 <- MCMCglmm(
  DomZooFormula1,
  data = Viruses,
  random =~ vFamily,
  prior = DomZooPrior2,
  family = "categorical",
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)

DomZooMultiv2 <- MCMCglmm(
  cbind(PropDomestic, IsZoonotic) ~ trait - 1 + trait:(Records),
  prior = DomZooMultivPrior2,
  rcov =~ us(trait):units,
  data = Viruses,
  random =~ us(trait):vFamily,
  family = c("gaussian","categorical"),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)

ZooDomSel <- INLAModelSel("Human", Viruscovar, "vFamily", "iid", "binomial", Viruses)

Attempt <- inla(as.formula(paste0("Human ~ ", paste(ZooDomSel$Removed[[6]], collapse = " + "), " + f(vFamily, model = 'iid') + f(PropDomestic, model = 'rw1')")),
                data = Viruses, 
                family = "binomial", 
                control.compute = list(dic = T))
