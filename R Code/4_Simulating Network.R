
# Generating the network of viruses/hosts

load("ZI_runs.Rdata")

i = 1

N = dim(FinalHostMatrix)[1]

CountColumns <- list(1:7*2-1, 15:662)
ZIColumns <- list(1:7*2, 663:1310)

RowSampled <- sample(1:nrow(ModelList[[i]]$Sol), 1)

CountFXSample <- ModelList[[i]]$Sol[RowSampled, unlist(CountColumns)]
ZIFXSample <- ModelList[[i]]$Sol[RowSampled, unlist(ZIColumns)]

XZMatrix <- cbind(ModelList[[i]]$X, ModelList[[i]]$Z)

CountXZMatrix <- as.matrix(XZMatrix[,unlist(CountColumns)])
ZIXZMatrix <- as.matrix(XZMatrix[,unlist(ZIColumns)])

CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))

CountOutput <- CountOutput[1:N]
ZIOutput <- ZIOutput[(N+1):(2*N)]

Responses <- cbind(ZIOutput, CountOutput)

# Trying it without random effects ####

i = 2

PredList <- list()

for(x in 1:1000){
  
  RowSampled <- sample(1:nrow(ModelList[[i]]$Sol), 1)
  
  CountFXSample <- ModelList[[i]]$Sol[RowSampled, CountColumns[[1]]]
  ZIFXSample <- ModelList[[i]]$Sol[RowSampled, ZIColumns[[1]]]
  
  XZMatrix <- cbind(ModelList[[i]]$X, ModelList[[i]]$Z)
  
  CountXZMatrix <- as.matrix(XZMatrix[,CountColumns[[1]]])
  ZIXZMatrix <- as.matrix(XZMatrix[,ZIColumns[[1]]])
  
  CountOutput <- c(CountFXSample %*% t(CountXZMatrix))
  ZIOutput <- c(ZIFXSample %*% t(ZIXZMatrix))
  
  CountOutput <- CountOutput[1:N]
  ZIOutput <- ZIOutput[(N+1):(2*N)]
  
  Responses <- cbind(ZIOutput, CountOutput)
  
  PZero <- logit(ZIOutput)
  PCount <- exp(CountOutput)*(1-PZero)
  
  PredList[[x]] <- PCount
  
}
