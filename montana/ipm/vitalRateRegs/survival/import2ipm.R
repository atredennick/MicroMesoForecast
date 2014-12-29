# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Test script to pull in growth parameters for IPM #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#Several arguments are required
#' @param doYear  Specific climate year to pull year random effects
#' @param mcDraw  A numeric scalar or vector for the row(s) of MCMC to draw parameters from 
#' @param doSpp   A character scalar for the current speicies
# @param group   A numeric scalar with group indicator

library(reshape2)
library(plyr)

sppList=c("ARTR","HECO","POSE","PSSP") #in here for now; remove after testing

####
#### Read in full MCMC output and format as data frame
####
MCMC <- readRDS("survivalParamsMCMC.rds")
psurv2 <- melt(MCMC)
psurv2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=5))
psurv2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*13),
                 rep("nb", times=4*3000),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(psurv2)[1] <- "Iter"
psurv <- psurv2[,c(1,3:5)]; rm(psurv2)
psurvAll <- subset(psurv, Coef=="gInt"|Coef=="nb"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2")
psurvYrs <- subset(psurv, Coef=="beta" | Coef=="intYr")
years <- unique(allD$year)[2:14]+1900
psurvYrs$Year <- c(rep(rep(years, each=3000), each=4),
                   rep(rep(years, each=3000), each=4))

####
#### Now get subset defined by parameters; in a function
####
getSurvCoefs <- function(doYear, mcDraw, doSpp){
  #First do coefficients with random year effects
  survNowYr <- subset(psurvYrs, Spp==doSpp & Year==doYear & Iter==mcDraw)
  iID <- which(survNowYr$Coef=="intYr")
  intercept <- survNowYr$value[iID]
  sID <- which(survNowYr$Coef=="beta")
  size <- survNowYr$value[sID]
  
  #Now do climate and competition fixed effects
  survNow <- subset(psurvAll, Spp==doSpp & Iter==mcDraw)
  cID <- which(survNow$Coef=="rain1"|survNow$Coef=="rain2"|survNow$Coef=="temp1"|survNow$Coef=="temp2")
  climEffs <- survNow$value[cID]
  dd <- survNow$value[which(survNow$Coef=="nb")]
  
  #Collate all parameters for output
  Spars=list(intcpt=intercept, 
             slope=size,
             nb=dd,
             clim=climEffs)
  
  out(Spars)
}
