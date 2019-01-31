
# Trivariate model ####

library(MCMCglmm)

NARows <-function(df, vars){
  apply(df[,vars], 1, function(a){
    any(is.na(a)|a=="Inf"|a=="-Inf")
  })
}

Viruscovar <- c("vDNAoRNA","vEnvelope","vVectorYNna","vPubMedCites","HostRangeMean")[1]
Viruscovar <- c("HostRangeMean")[1]

TestViruses <- Viruses[-which(NARows(Viruses, Viruscovar)),]
TestViruses <- Viruses[-which(is.na(Viruses$HostRangeMean)),]

TrivPrior1 <- 
  list(R = list(V = diag(3), nu = 0.002, fix = 1))

TrivPrior1g <- 
  list(R = list(V = diag(3), nu = 0.002))

TrivPrior2 <- 
  list(R = list(V = diag(3), nu = 0.002, fix = 1),
       G = list(G1 = diag(3), nu = 0.002))

TrivFormula <- as.formula(paste0("cbind(Human, Domestic, Wildlife)", " ~ - 1 + trait + trait:(", paste(c(Viruscovar), collapse = " + "), ")"))

mf = 1

TrivMCMC1 <- MCMCglmm(
  TrivFormula,
  prior = TrivPrior1,
  rcov =~ us(trait):units,
  data = TestViruses,
  #random =~ us(trait):vFamily,
  family = rep("categorical", 3),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)

TrivMCMC2 <- MCMCglmm(
  TrivFormula,
  prior = TrivPrior2,
  rcov =~ us(trait):units,
  data = TestViruses,
  random =~ us(trait):vFamily,
  family = rep("categorical", 3),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)

TrivMCMC1g <- MCMCglmm(
  TrivFormula,
  prior = TrivPrior1g,
  rcov =~ us(trait):units,
  data = TestViruses,
  #random =~ us(trait):vFamily,
  family = rep("gaussian", 3),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)

TrivMCMC2g <- MCMCglmm(
  TrivFormula,
  prior = TrivPrior2,
  rcov =~ us(trait):units,
  data = TestViruses,
  random =~ us(trait):vFamily,
  family = rep("gaussian", 3),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)

# Trying it with counts ####

TrivFormula2 <- as.formula(paste0("cbind(DomesticCount, WildlifeCount, Human)", " ~ - 1 + trait + (", paste(c(Viruscovar), collapse = " + "), "):trait"))

TrivPrior3 <- 
  list(R = list(V = diag(3), nu = 0.002, fix = 3))

TrivMCMCCount <- MCMCglmm( # This doesn't make sense
  TrivFormula2,
  prior = TrivPrior3,
  rcov =~ us(trait):units,
  data = TestViruses,
  #random =~ us(trait):vFamily,
  family = c("poisson", "poisson", "categorical"),
  nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
  thin = 10*mf, burnin = 3000*mf
)
