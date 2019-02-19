
library('mgcv')
library('brms')
library('ggplot2')
library('schoenberg')
library('bayesplot')

DNAGAM <- readRDS("DNAGAM.rds")


schoenberg::draw(DNAGAM)

plot(DNAGAM)

pp_check(DNAGAM, nsamples = 30)
pp_check(DNAGAM, nsamples = 30, type = "ecdf_overlay")


