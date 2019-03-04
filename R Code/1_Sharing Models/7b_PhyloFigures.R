
# Importing Phylopics ####

library(rphylopic); library(tidyverse)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)

lapply(levels(Panth1$hOrder), function(a){
  
  Species <- Panth1 %>% filter(hOrder==a) %>% select(Sp) 
  
  Species$Sp <- Species$Sp %>% str_replace_all("_" , " ")
  
  lapply(Species$Sp, function(b){
    
    pic <- rphylopic::search_text(text = b, options = "names")
    
    output <- search_images(pic[1], options=c("pngFiles", "credit", "canonicalName"))
    
    pics <- image_data(output[4,], size = "256")
    
    save_png(pics, target = paste0("Figures/",#a,"/",
                                   b,".png"))
    
  })
  
})



person <- name_search(text = "Homo sapiens", options = "namebankID")
pig <- name_search(text = "Sus scrofa", options = "namebankID")

img <- image_data("27356f15-3cf8-47e8-ab41-71c6260b2724", size = "512")

x <- search_text(text = "Homo sapiens", options = "names")
output <- search_images(x[1], options=c("pngFiles", "credit", "canonicalName"))
pic = image_data(output, size = "256")

load("Output Files/Panth1.Rdata")

p = BarGraph(Panth1, "hOrder", "AllPredDegree") + 
  add_phylopic(PNGList[[1]], 1)


# Doing this manually

library(png); library(grid)

PNGFiles <- list.files("Phylopics")
PNGList <- lapply(PNGFiles, function(a) readPNG(source = paste0("Phylopics/",a)))
PNGNames <- PNGFiles %>% substr(1,nchar(.)-4)
names(PNGList) <- PNGNames

g = rastergrob(PNGList[[1]])

p = BarGraph(Panth1, "hOrder", "AllPredDegree") 

p + rasterImage(PNGList[[1]],
                xleft = which(df2$hOrder)==2)

df = Panth1
x = "hOrder" 
y = "AllPredDegree" 
z = x
geom = "Bar" 
text = "N"
labels = NA
order = T 
Just = T

require(ggplot2)
require(reshape2)

df2 <- data_summary(droplevels(na.omit(df[, c(x, y, z)])), 
                    varname = y, groupnames = unique(c(x, z)))
df2 <- df2[order(df2[, y]), ]
df2[, x] <- factor(df2[, x], levels = unique(df2[, x]))

if (text %in% c("Prevalence", "N")) {
  p <- ggplot(df2, aes(as.numeric(df2[, x]), df2[, y], fill = as.factor(df2[, 
                                                                            z]))) + geom_bar(stat = "identity", colour = "black", 
                                                                                             position = position_dodge(0.9)) + geom_errorbar(aes(ymin = df2[, 
                                                                                                                                                            y] - se, ymax = df2[, y] + se), width = 0.2, 
                                                                                                                                             position = position_dodge(0.9)) + THEME + xlab(x) + 
    ylab(y) + labs(fill = z) + geom_text(aes(y = ifelse(df2[, 
                                                            y] + se > 0, df2[, y] + se + PositionT, PositionT), 
                                             label = N), position = position_dodge(0.9)) + 
    geom_text(aes(label = labels, y = ifelse(df2[, 
                                                 y] + se > 0, df2[, y] + se + PositionT * 1.5, 
                                             PositionT * 1.5)), position = position_dodge(0.9)) + 
    theme(axis.text.x = element_text(angle = Angle, 
                                     hjust = Hjust))
}

for(q in PNGNames){
  
  g = rasterGrob(PNGList[[q]])
  
  OrderWhich = which(df2$hOrder == q)
  
  p = p + annotation_custom(g,
                            xmin = OrderWhich - 0.5,
                            xmax = OrderWhich + 0.5,
                            ymin = 5,
                            ymax = 15) #df2[OrderWhich,"AllPredDegree"])
  
}

p + scale_x_discrete(labels = levels(df2$hOrder),
                       limits = 1:nlevels(df2$hOrder))

p + scale_x_discrete(labels = levels(df2$hOrder),
                     limits = 1:nlevels(df2$hOrder))+ 
  ggsave("Figures/Phylopic.jpeg", units = "mm", height = 120, width = 200)
