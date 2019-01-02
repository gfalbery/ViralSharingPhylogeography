
# Centrality Model Selection: Viruses

library(INLA)

Viruscovar <- c("Records","vCytoReplicTF","vDNAoRNA","vEnvelope","vGenomeAveLengthLn","vPubMedCitesLn","vSegmentedTF","vSSoDS","vVectorYNna")
Infects <- c("Wildlife", "Domestic", "Human")

Hostcovar <- c("UrbRurPopRatioLn","AreaHost","HabInhabitedChgLn","hAllZACitesLn","hMarOTerr","hMassGramsPVR","HumPopDensLnChg","TotHumPopLn")
Domcovar <- c("domestic_category","DOMYearBP")

ViralCentralitySel <- list()

ViralCentralitySel[[1]] <- INLAModelSel(vCentrality[1], Viruscovar, "vFamily", "iid", vDists[1], Viruses)
ViralCentralitySel[[2]] <- INLAModelSel(vCentrality[2], Viruscovar, "vFamily", "iid", vDists[2], Viruses)
ViralCentralitySel[[3]] <- INLAModelSel(vCentrality[3], Viruscovar, "vFamily", "iid", vDists[3], Viruses)

ViralCentralityAdd <- list()

ViralCentralityAdd[[1]] <- INLAModelAdd(vCentrality[1], 1, c(Viruscovar, Infects), "vFamily", "iid", vDists[1], Viruses)
