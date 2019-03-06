
library(mgcv); library(MASS)


for(i in c("Space", "Phylo")){
  
  lp <- predict(BAMList[[Resps[r]]], newdata = FitList[[Resps[r]]] %>% 
                  filter(i == last(unique(i))), 
                type = "lpmatrix") %>% 
    as.data.frame()
  
  coefs <- coef(BAMList[[Resps[r]]])
  vc <- vcov(BAMList[[Resps[r]]])
  
  sim <- mvrnorm(100, mu = coefs, Sigma = vc)
  
  want <- lp %>% colnames
  
  lp <- lp %>% as.matrix #%>% logistic
  
  fits <- lp[, want] %*% t(sim[, want]) %>% as.data.frame() %>%
    mutate(i = FitList[[Resps[r]]][,i])
  
  PostList[[Resps[r]]][[i]] <- gather(fits, key = "Draw", value = "Fit", -Space)
  
}

ggplot(FitLong, aes(Space, logistic(Fit), colour = Draw)) + geom_line(alpha = 0.3)


