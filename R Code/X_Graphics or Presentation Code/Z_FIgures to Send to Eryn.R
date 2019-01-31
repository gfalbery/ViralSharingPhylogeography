
# Figs for Eryn #####

ModelNames2 <- ModelNames[c(1,3,4,7)]
Resps2 <- Resps[1:3]

jpeg("Figures/Space and Phylogeny Matrices Improve Model Fit.jpeg", units = "mm", width = 200, height = 200, res = 300)
lapply(1:length(CentralityList[Resps2]), function(a){ 
  INLADICFig(CentralityList[[a]][ModelNames2], ModelNames = ModelNames2) + 
    ggtitle(Resps2[a]) + theme(legend.position = "none", axis.text.x = element_text(angle = 10, hjust = 1)) + labs(x = "Model")
}) %>% 
  arrange_ggplot2(nrow = 3)
dev.off()

jpeg("Figures/Records and Eigenvector Display Suspiciously Good Fit.jpeg", units = "mm", width = 200, height = 100, res = 300)

lapply(1:length(CentralityList[Resps2]), function(a){ 
  qplot(TestHosts[,Resps2[a]],
        CentralityList[[a]][[7]]$summary.fitted.values$mean[1:dim(TestHosts)[1]]) + 
    ggtitle(Resps2[a]) + labs(x = paste("Data", Resps2[a]), y = paste("Fitted", Resps2[a]))
}) %>%
  arrange_ggplot2(nrow = 1)

dev.off()
