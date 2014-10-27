# Get the DIC ranking tables and add the formula as a column

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

####
#### Growth
####
source("../vitalRateRegressions/fixedEffectsModels_Growth.R")
tabG <- read.csv("../vitalRateRegressions/topModels_Growth.csv")
index <- tabG$modelNum
modVec <- fixed.forms[index]
tabG$Model <- as.vector(modVec)


####
#### Survival
####
source("../vitalRateRegressions/fixedEffectsModels_Survival.R")
tabS <- read.csv("../vitalRateRegressions/topModels_Survival.csv")
index <- tabS$modelNum
modVec <- fixed.forms[index]
tabS$Model <- as.vector(modVec)


####
#### Colonization
####
source("../vitalRateRegressions/fixedEffectsModels_Colonization.R")
tabC <- read.csv("../vitalRateRegressions/topModels_Colonization.csv")
index <- tabC$modelNum
modVec <- fixed.forms[index]
tabC$Model <- as.vector(modVec)

