
# Calculating MCMC repeatability for a binomial model #####

MCBinRep <- function(Model, scale = "original"){
  
  mat <- matrix(NA, nrow = dim(Model$VCV)[2], ncol = 4)
  
  if(scale == "original"){
    
    Beta0 <- sapply(1:dim(Model$Sol)[1], function(z) mean(as.matrix(Model$X) %*% 
                                                            as.matrix(Model$Sol[z, 1:ncol(Model$X)])))
    
    P = exp(Beta0)/(1 + exp(Beta0))
    
    for (j in 1:dim(Model$VCV)[2]) {
      Va <- Model$VCV[, j]
      Ve <- rowSums(Model$VCV)
      Repeatability1 <- ((Va*P^2)/(1 + exp(Beta0))^2)/(((Ve)*P^2)/((1 + exp(Beta0))^2 + P*(1 - P)))
      
      mat[j, ] <- c(colnames(Model$VCV)[j], round(posterior.mode(Repeatability1), 
                                                  digits = 2), 
                    round(HPDinterval(Repeatability1)[1], digits = 2), 
                    round(HPDinterval(Repeatability1)[2],2))    
    }
  }
  
  if(scale == "link"){
    for (j in 1:dim(Model$VCV)[2]) {
      Va <- Model$VCV[, j]
      Ve <- rowSums(Model$VCV)
      Repeatability1 <- (Va)/(Ve + (pi^2/3))
      mat[j, ] <- c(colnames(Model$VCV)[j], round(posterior.mode(Repeatability1), 
                                                  digits = 2), 
                    round(HPDinterval(Repeatability1)[1], digits = 2), 
                    round(HPDinterval(Repeatability1)[2],2))
    } 
    
  }
  
  colnames(mat) <- c("Component", "Mode", "lHPD", "uHPD")
  data.frame(mat)
  
}




MCMCRep
function (Model, scale = "original") 
{
  require(MCMCglmm)
  if (scale == "original") {
    if (unique(Model$family) == "poisson") {
      Beta0 <- sapply(1:dim(Model$Sol)[1], function(z) mean(as.matrix(Model$X) %*% 
                                                              as.matrix(Model$Sol[z, 1:ncol(Model$X)])))
      mat <- matrix(NA, nrow = dim(Model$VCV)[2], ncol = 4)
      for (j in 1:dim(Model$VCV)[2]) {
        Va <- Model$VCV[, j]
        Ve <- rowSums(Model$VCV)
        Expected <- exp(Beta0 + (0.5 * (Ve)))
        Repeatability1 <- (Expected * (exp(Va) - 1))/(Expected * 
                                                        (exp(Ve) - 1) + 1)
        mat[j, ] <- c(colnames(Model$VCV)[j], round(posterior.mode(Repeatability1), 
                                                    digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                   1], digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                                                                      2], digits = 2))
      }
    }
    if (unique(Model$family) == "gaussian") {
      mat <- matrix(NA, nrow = dim(Model$VCV)[2], ncol = 4)
      for (j in 1:dim(Model$VCV)[2]) {
        Va <- Model$VCV[, j]
        Ve <- rowSums(Model$VCV)
        Repeatability1 <- (Va/Ve)
        mat[j, ] <- c(colnames(Model$VCV)[j], round(posterior.mode(Repeatability1), 
                                                    digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                   1], digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                                                                      2], digits = 2))
      }
    }
  }
  if (scale == "link") {
    if (unique(Model$family) == "poisson") {
      Beta0 <- sapply(1:dim(Model$Sol)[1], function(z) mean(as.matrix(Model$X) %*% 
                                                              as.matrix(Model$Sol[z, 1:ncol(Model$X)])))
      mat <- matrix(NA, nrow = dim(Model$VCV)[2], ncol = 4)
      for (j in 1:dim(Model$VCV)[2]) {
        Va <- Model$VCV[, j]
        Ve <- rowSums(Model$VCV)
        Repeatability1 <- Va/(Ve + log(1/exp(Beta0) + 
                                         1))
        mat[j, ] <- c(colnames(Model$VCV)[j], round(posterior.mode(Repeatability1), 
                                                    digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                   1], digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                                                                      2], digits = 2))
      }
    }
    if (unique(Model$family) == "gaussian") {
      mat <- matrix(NA, nrow = dim(Model$VCV)[2], ncol = 4)
      for (j in 1:dim(Model$VCV)[2]) {
        Va <- Model$VCV[, j]
        Ve <- rowSums(Model$VCV)
        Repeatability1 <- (Va/Ve)
        mat[j, ] <- c(colnames(Model$VCV)[j], round(posterior.mode(Repeatability1), 
                                                    digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                   1], digits = 2), round(HPDinterval(Repeatability1)[, 
                                                                                                                                                      2], digits = 2))
      }
    }
  }
  colnames(mat) <- c("Component", "Mode", "lHPD", "uHPD")
  data.frame(mat)
}
