
# 0i_Importing Ecological data ####

EltonTraits <- read.delim("data/MamFuncDat.txt")
EltonTraits$Scientific <- EltonTraits$Scientific %>% str_replace(" ","_")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt")
Panth2 <- read.delim("data/PanTHERIA_1-0_WR93_Aug2008.txt")

Panth1$MSW05_Binomial <- 
  Panth1$MSW05_Binomial %>% str_replace(" ", "_")

EcoDistVars <- c(
  "X5.1_AdultBodyMass_g", 
  "X15.1_LitterSize", 
  "X16.1_LittersPerYear",
  "X22.1_HomeRange_km2",
  "X27.2_HuPopDen_Mean_n.km2",
  "X28.1_Precip_Mean_mm",
  "X28.2_Temp_Mean_01degC",
  "X6.2_TrophicLevel",
  "X10.2_SocialGrpSize",
  "X6.1_DietBreadth", 
  "X12.1_HabitatBreadth"
)

kt <- ktab.data.frame(Panth1[,EcoDistVars])
