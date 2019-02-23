
# Doing this to validate rather than predict ####

# Rscript "R Code/1_Sharing Models/4_Host Validation.R"

library(tidyverse); library(parallel); library(ggregplot); library(ape)

source("R Code/00_Master Code.R")
load("AllSimsGAM.Rdata")

if(file.exists("Output Files/GAMValidation.Rdata")) load("Output Files/GAMValidation.Rdata")

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

GAMValidation <- list()

print("Start Validating!")

a = 1

for(a in a:length(VirusAssocs)){
  
  print(names(VirusAssocs)[a])
  
  pHosts <- VirusAssocs[[a]]
  
  pHosts <- intersect(pHosts, AllMammals)
  
  if(length(pHosts)>0){
    
    EstList <- parallel::mclapply(1:length(AllSims), function(x){
      
      FocalNet <- AllSims[[x]] %>% as.matrix
      
      pHosts2 <- intersect(pHosts, rownames(FocalNet))
      
      ValidEst <- list()
      
      for(b in pHosts2){
        
        pHosts4 <- setdiff(pHosts2, b)
        
        pHosts3 <- setdiff(colnames(FocalNet), pHosts4)
        
        Estimates <- FocalNet[pHosts4, pHosts3]
        
        if(is.null(dim(Estimates))) Estimates <- rbind(Estimates, Estimates)/2
        
        Ests <- data.frame(Sp = names(sort(colSums(Estimates))),
                           Count = sort(colSums(Estimates))/nrow(Estimates),
                           Iteration = x) %>%
          mutate(Focal = ifelse(Sp==b, 1, 0))
        
        rownames(Ests) <- Ests$Sp
        
        ValidEst[[b]] <- Ests
        
      }
      
      ValidEst
      
    }, mc.cores = 40)
    
    GAMValidation[[names(VirusAssocs)[a]]] <- EstList
    
  } else GAMValidation[[names(VirusAssocs)[a]]] <- NA
  
  if(a %% 10 == 0){
    
    GAMValid <- GAMValidation %>% lapply(function(a){
      
      if(!is.null(names(a[[1]]))){
        
        b = map(names(a[[1]]), function(b) map(a, b) %>% bind_rows) %>% bind_rows
        
        c = b %>% group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>% 
          slice(order(Count, decreasing = T)) %>%
          mutate(Focal = factor(Focal))
        
      } else c = NA
      
      return(c)
    })
    
    save(GAMValid, file = "Output Files/GAMValidation.Rdata")
    
  }
}

GAMValid <- GAMValidation %>% lapply(function(a){
  
  if(!is.null(names(a[[1]]))){
    
    b = map(names(a[[1]]), function(b) map(a, b) %>% bind_rows) %>% bind_rows
    
    c = b %>% group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>% 
      slice(order(Count, decreasing = T)) %>%
      mutate(Focal = factor(Focal))
    
  } else c = NA
  
  return(c)
})

save(GAMValid, file = "Output Files/GAMValidation.Rdata")

load("Output Files/GAMValidation.Rdata")
load("Output Files/ModelValidation.Rdata")

KeepPredictions <- (1:length(GAMValid))[-which(sapply(GAMValid, function(a) any(is.na(a))))]

FocalRank <- function(x){
  
  y <- x[x[,"Focal"]==1,"Count"]
  z <- x[x[,"Focal"]==0,"Count"]
  
  (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
  
}

GAMValidSummary <- data.frame(
  
  Virus = names(GAMValid)[KeepPredictions],
  
  NHosts = map(GAMValid[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
  
  No = KeepPredictions,
  
  MeanRank = sapply(GAMValid[KeepPredictions], function(a) mean(FocalRank(a)))
  
) %>% slice(order(MeanRank))


ValidSummary <- data.frame(
  
  Virus = names(Valid)[KeepPredictions],
  
  NHosts = map(Valid[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
  
  No = KeepPredictions,
  
  MeanRank = sapply(Valid[KeepPredictions], function(a) mean(FocalRank(a)))
  
) #%>% slice(order(MeanRank))


CompSummary <- merge(GAMValidSummary, ValidSummary[,c("Virus","MeanRank")], by = "Virus", suffixes = c(".GAM",".GLM"), all.x = T)

CompSummary %>% ggplot(aes(MeanRank.GAM, MeanRank.GLM)) + geom_point() + geom_smooth() + geom_abline() + coord_fixed()


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

rownames(GAMValidSummary) <- GAMValidSummary$Virus

Viruses[,"PredictionSuccess"] <- GAMValidSummary[as.character(Viruses$Sp), "MeanRank"]

PredHostPlot <- function(virus, threshold = 10, focal = c(1,0), facet = FALSE){
  
  require(ggtree)
  
  Df <- GAMValid[[virus]] %>% mutate(Focal = as.numeric(as.character(Focal)))
  Df$Include <- ifelse(Df$Focal%in%focal, 1, 0)
  Df <- Df %>% filter(Include==1)
  
  Df$Rank = nrow(Df) - rank(Df$Count)
  
  if(all(focal == 1)) threshold <- nrow(Df)
  
  PredHosts <- Df %>% filter(Rank < threshold) #%>% select(Sp) %>% unlist
  if(length(focal)==2) PredHosts <- Df %>% filter(Rank < threshold|Focal==1) #%>% select(Sp) %>% unlist
  
  PredHostPolygons <- FullPolygons %>% filter(Host%in%PredHosts$Sp) %>% 
    left_join(PredHosts, by = c("Host" = "Sp")) %>%
    mutate(Host = factor(Host, levels = Df[order(Df$Rank, Df$Focal, decreasing = TRUE),"Sp"] %>% unlist))
  
  VirusName <- str_replace_all(virus, "_", " ")
  
  Focal <- GAMValid[[virus]] %>% mutate(Focal = as.numeric(as.character(Focal))) %>% filter(Focal==1) %>% select(Sp) %>% unlist
  
  if(0 %in%focal) Predicted <- PredHosts %>% filter(Focal==0) %>% select(Sp) %>% unlist else{
    Predicted <- NULL
  }
  
  Groups <- ifelse(STFull$tip.label%in%Focal,"Known",ifelse(STFull$tip.label%in%Predicted,"Predicted",""))
  
  groupInfo <- split(STFull$tip.label, Groups)
  chiroptera <- groupOTU(STFull, groupInfo)
  
  plot2 <- ggtree(chiroptera, aes(color=group, alpha = group)) +
    scale_colour_manual(values = c("black", "red", "blue")) +
    scale_alpha_manual(values = c(0.01,0.5,1)) +
    theme(legend.position = "none")
  
  g2 = ggplotGrob(plot2)
  
  xmin = -1.9*(10^7) 
  xmax = -1.2*(10^7)
  ymin = -9*(10^6)
  ymax = 9*(10^6)
  
  rectborder <- 60000
  rect = data.frame(long = c(xmin-rectborder, xmin-rectborder, xmax+rectborder, xmax+rectborder), 
                    lat = c(ymin-rectborder, ymax+rectborder, ymax+rectborder, ymin-rectborder))
  
  if(facet == FALSE){
    
    ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) + 
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) + # alpha = max(Rank)-Rank)) +
      labs(#alpha = "Inverse Rank", 
        title = paste(ifelse(length(focal)==2, "All", ifelse(focal==1, "Known", "Predicted")), VirusName, "Hosts")) +
      coord_fixed() +
      #theme(legend.position = "none") +
      scale_x_continuous(breaks = -10:10*2000000) +
      scale_y_continuous(breaks = -5:5*2000000) +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat), 
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2, 
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax) %>% return #facet_wrap(~Host)
    
  }else{
    
    ggplot(PredHostPolygons, aes(long, lat, group = paste(Host, group))) + 
      geom_polygon(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), fill = "white", colour = "black") +
      geom_polygon(fill = NA, aes(colour = Host)) + # alpha = max(Rank)-Rank)) +
      labs(#alpha = "Inverse Rank", 
        title = paste(ifelse(length(focal)==2, "All", ifelse(focal==1, "Known", "Predicted")), VirusName, "Hosts")) +
      coord_fixed() +
      #theme(legend.position = "none") +
      scale_x_continuous(breaks = -10:10*2000000) +
      scale_y_continuous(breaks = -5:5*2000000) +
      geom_polygon(inherit.aes = F, data = rect, aes(long, lat), 
                   fill = "grey", colour = "grey") +
      annotation_custom(grob = g2, 
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax) + 
      facet_wrap(~Focal, ncol = 1, 
                 labeller = labeller(Focal = c("0" = "Predicted", "1" ="Known"))) %>% return #facet_wrap(~Host)
  }
}

pdf("HostPredictions.pdf", width = 14, height = 6)

lapply(Viruses$Sp[KeepPredictions], function(a) PredHostPlot(a, focal = 0)) %>% return

dev.off()

pdf("HostKnown.pdf", width = 14, height = 6)

lapply(Viruses$Sp[KeepPredictions], function(a) PredHostPlot(a, focal = 1)) %>% return

dev.off()

pdf("HostBoth.pdf", width = 14, height = 12)

lapply(Viruses$Sp[KeepPredictions], function(a) PredHostPlot(a, focal = c(0,1), facet = TRUE)) %>% return

dev.off()
