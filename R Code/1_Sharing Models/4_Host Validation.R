
# Doing this to validate rather than predict ####

library(tidyverse); library(parallel); library(ggregplot)

source("R Code/00_Master Code.R")
load("AllSims.Rdata")

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]

HostValidation <- list()

print("Start Validating!")

a = length(Valid)

for(a in a:(length(VirusAssocs))){
  
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
    
    HostValidation[[names(VirusAssocs)[a]]] <- EstList
    
  } else HostValidation[[names(VirusAssocs)[a]]] <- NA
  
  if(a %% 10 == 0){
    
    Valid <- HostValidation %>% lapply(function(a){
      
      if(!is.null(names(a[[1]]))){
        
        b = map(names(a[[1]]), function(b) map(a, b) %>% bind_rows) %>% bind_rows
        
        c = b %>% group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>% 
          slice(order(Count, decreasing = T)) %>%
          mutate(Focal = factor(Focal))
        
      } else c = NA
      
      return(c)
    })
    
    #save(Valid, file = "ModelValidation.Rdata")
    
  }
}

# Trying it out ####

Valid <- HostValidation %>% lapply(function(a){
  
  if(!is.null(names(a[[1]]))){
    
    b = map(names(a[[1]]), function(b) map(a, b) %>% bind_rows) %>% bind_rows
    
    c = b %>% group_by(Sp, Focal) %>% dplyr::summarise(Count = mean(Count)) %>% slice(order(Count, decreasing = T)) %>%
      mutate(Focal = factor(Focal))
    
  } else c = NA
  
  return(c)
})

#save(Valid, file = "ModelValidation.Rdata")

KeepPredictions %>% 
  lapply(function(a) ggplot(Valid[[a]], aes(Focal, Count, colour = Focal)) + 
           ggforce::geom_sina() + theme(legend.position = "none") +
           ggtitle(names(VirusAssocs)[[a]])) %>% 
  arrange_ggplot2(ncol = 5)

FocalRank <- function(x){
  
  y <- x[x[,"Focal"]==1,"Count"]
  z <- x[x[,"Focal"]==0,"Count"]
  
  (length(z$Count) + 2) - sapply(y$Count, function(a) rank(c(a,z$Count))[1])
  
}

load("ModelValidation.Rdata")

KeepPredictions <- (1:length(Valid))[-which(sapply(Valid, function(a) any(is.na(a))))]

ValidSummary <- data.frame(
  
  Virus = names(VirusAssocs)[KeepPredictions],
  
  NHosts = map(Valid[KeepPredictions], "Focal") %>% sapply(function(a) which(a=="1") %>% length),
  
  No = KeepPredictions,
  
  MeanRank = sapply(Valid[KeepPredictions], function(a) mean(FocalRank(a)))
  
) %>% slice(order(MeanRank))

rownames(ValidSummary) <- ValidSummary$Virus

Viruses[,"PredictionSuccess"] <- ValidSummary[as.character(Viruses$Sp), "MeanRank"]

ggplot(ValidSummary, aes(log10(NHosts), log10(MeanRank))) + geom_smooth() + geom_text(aes(label = Virus))

qplot(log10(ValidSummary$MeanRank))

table(cut(ValidSummary$MeanRank, breaks = (10^(0:4))-1, labels = c("<10", "<100", "<1000", ">1000")))

ggplot(Viruses, aes(HostRangeMean, log(PredictionSuccess+1))) + geom_text(aes(label = Sp)) + geom_smooth()
ggplot(Viruses, aes(HostRangeMax, log(PredictionSuccess+1))) + geom_text(aes(label = Sp)) + geom_smooth()
ggplot(Viruses, aes(HostRangeMin, log(PredictionSuccess+1))) + geom_text(aes(label = Sp)) + geom_smooth()

PredHostPlot <- function(virus, threshold = 10, focal = c(1,0), facet = FALSE){
  
  Df <- Valid[[virus]] %>% mutate(Focal = as.numeric(as.character(Focal)))
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
  
  Focal <- Valid[[virus]] %>% mutate(Focal = as.numeric(as.character(Focal))) %>% filter(Focal==1) %>% select(Sp) %>% unlist
  
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
      geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), colour = "black") +
      geom_polygon(aes(fill = Host, alpha = max(Rank)-Rank)) +
      labs(alpha = "Inverse Rank", 
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
      geom_path(data = WorldMap, inherit.aes = F, aes(long, lat, group = group), colour = "black") +
      geom_polygon(aes(fill = Host, alpha = max(Rank)-Rank)) +
      labs(alpha = "Inverse Rank", 
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

