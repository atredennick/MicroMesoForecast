####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#### Vital Rate Regression MCMC wrapper ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

####
#### Growth
####
setwd("./growth/")
source("growthWrapper.R")

####
#### Survival
####
setwd("../survival/")
source("survivalWrapper.R")

####
#### Colonization
####
setwd("../recruitment/")
source("recruitmentWrapper.R")
