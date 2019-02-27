# 4b_Plotting Host Validation Results ###

# Plotting results ####

KeepPredictions %>% 
  lapply(function(a) ggplot(GAMValid[[a]], aes(Focal, Count, colour = Focal)) + 
           ggforce::geom_sina() + theme(legend.position = "none") +
           ggtitle(names(GAMValid)[[a]])) %>% 
  arrange_ggplot2(ncol = 3)

KeepPredictions %>% 
  lapply(function(a) ggplot(Valid[[a]], aes(Focal, Count, colour = Focal)) + 
           ggforce::geom_sina() + theme(legend.position = "none") +
           ggtitle(names(Valid)[[a]])) %>% 
  arrange_ggplot2(ncol = 3)

# PDF

Viruses[,"PredictionSuccess"] <- GAMValidSummary[as.character(Viruses$Sp), "MeanRank"]

pdf("HostPredictions.pdf", width = 14, height = 6)

lapply(Viruses$Sp[KeepPredictions], function(a) PredHostPlot(a, focal = 0)) %>% return

dev.off()

pdf("HostKnown.pdf", width = 14, height = 6)

lapply(Viruses$Sp[KeepPredictions], function(a) PredHostPlot(a, focal = 1)) %>% return

dev.off()

pdf("HostBoth.pdf", width = 14, height = 12)

lapply(Viruses$Sp[KeepPredictions], function(a) PredHostPlot(a, focal = c(0,1), facet = TRUE)) %>% return

dev.off()