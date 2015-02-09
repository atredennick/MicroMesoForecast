####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#### Vital Rate Regression MCMC wrapper ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

####
#### Growth
####
print("STARTING GROWTH")
setwd("./growth/")
# source("growthWrapper.R")
source("growthAllSpp_MCMC.R")
print("DONE WITH GROWTH")

####
#### Survival
####
print("STARTING SURVIVAL")
setwd("../survival/")
source("survivalAllSpp_MCMC.R")
print("DONE WITH SURVIVAL")

####
#### Colonization
####
print("STARTING RECRUITMENT")
setwd("../recruitment/")
source("recruitmentAllSpp_MCMC.R")
print("DONE WITH RECRUITMENT")
