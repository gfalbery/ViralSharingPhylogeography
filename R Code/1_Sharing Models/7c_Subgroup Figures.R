
Translate <- c("RNA Viruses",
               "DNA Viruses", 
               "Vector-Borne RNA Viruses",
               "Other RNA Viruses")

names(Translate) <- Resps[2:5]

Resps[2:5] %>% lapply(function(a){
  
  FitList[[a]] %>% 
    filter(Space == last(unique(Space))) %>%
    filter(DietSim == last(unique(DietSim))) %>%
    ggplot(aes(Phylo, Fit)) + geom_line() +
    lims(y = c(0,1)) +
    labs(y = "Predicted Viral Sharing", 
         x = "Phylogenetic Similarity", 
         title = a)
  
}) %>% arrange_ggplot2(ncol = 2)

jpeg("SIFigures/SubGroup_Space.jpeg", units = "mm", width = 150, height = 150, res = 300)

Resps[2:5] %>% lapply(function(a){
  
  PostList[[a]]$Space %>% 
    ggplot(aes(i, Fit, colour = Draw)) + geom_line(alpha = 0.3) + theme(legend.position = "none") +
    labs(x = "Space", y = "Model Estimate", title = Translate[a]) +
    geom_rug(data = DataList[[a]], inherit.aes = F, aes(x = Space), alpha = 0.01) +
    scale_colour_discrete_sequential(palette = AlberPalettes[2])
  
}) %>% arrange_ggplot2(ncol = 2)

dev.off()

jpeg("SIFigures/SubGroup_Phylo.jpeg", units = "mm", width = 150, height = 150, res = 300)

Resps[2:5] %>% lapply(function(a){
  
  PostList[[a]]$Phylo %>% 
    ggplot(aes(i, Fit, colour = Draw)) + geom_line(alpha = 0.3) + theme(legend.position = "none") +
    labs(x = "Phylo", y = "Model Estimate", title = Translate[a]) +
    geom_rug(data = DataList[[a]], inherit.aes = F, aes(x = Phylo), alpha = 0.01) +
    scale_colour_discrete_sequential(palette = AlberPalettes[1])
  
}) %>% arrange_ggplot2(ncol = 2)

dev.off()
