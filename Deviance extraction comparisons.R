
# Comparison of deviance models ####

load("~/Albersnet/Output Files/DevList.Rdata")
load("~/Albersnet/Output Files/TensorDevList.Rdata")

sapply(DevList, deviance)
sapply(TensorDevList, deviance)

OrigDev = deviance(DevList$FullModel)
RemoveDev = sapply(DevList[2:length(DevList)], deviance)

DevExplained = (RemoveDev - OrigDev)

(DevExplained/sum(DevExplained)) %>% round(2)

OrigDev = deviance(TensorDevList$FullModel)
RemoveDev = sapply(TensorDevList[2:length(TensorDevList)], deviance)

DevExplained = (RemoveDev - OrigDev)

(DevExplained/sum(DevExplained))# %>% round(2)
