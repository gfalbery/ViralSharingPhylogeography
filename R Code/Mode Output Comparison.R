
# Comparison of effects across different STAN models ####

RNABinModel <- readRDS("~/Albersnet/RNABinModel.rds")
DNABinModel <- readRDS("~/Albersnet/DNABinModel.rds")
VectorBinModel <- readRDS("~/Albersnet/VectorBinModel.rds")
NVectorBinModel <- readRDS("~/Albersnet/NVectorBinModel.rds")

list(
  stan_plot(RNABinModel, XBetas[2:4]),
  stan_plot(DNABinModel, XBetas[2:4]),
  stan_plot(VectorBinModel, XBetas[2:4]),
  stan_plot(NVectorBinModel, XBetas[2:4])
) %>% arrange_ggplot2

q$df[,XBetas]
