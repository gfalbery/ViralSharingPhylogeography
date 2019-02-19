# 0h_ Getting GLW in there ####

# Importing the AW files from Harvard Dataverse ####

GLWURL <- list(
  Buffalo = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/5U8MWI/T7QSGJ",
  Cattle = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/GIVQ75/I5CUJS",
  Chicken = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/SUFASB/AP1LHN",
  Duck = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/ICHCBH/NFEQ8V",
  Horse = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/7Q52MV/UAWH3Z",
  Goat = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/OCPH42/ZMVOOW",
  Pig = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/33N0JG/9LGGBS",
  Sheep = "https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/BLWPZN/XDIRM4"
)

LivestockNames <- names(GLWURL)

lapply(GLWURL, browseURL) # For some reason this is the only way this works 

GLWFiles <- paste0("data/GLW/",list.files("data/GLW"))
GLWRasters <- lapply(GLWFiles, raster::raster)

names(GLWRasters) <- LivestockNames
