# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Test script to pull in rect parameters for IPM #
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
MCMC <- readRDS("recruitmentParamsMCMC.rds")
prec2 <- melt(MCMC)
prec2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=5))
prec2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*13),
                 rep("nb", times=4*3000),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(prec2)[1] <- "Iter"
prec <- prec2[,c(1,3:5)]; rm(prec2)
precAll <- subset(prec, Coef=="gInt"|Coef=="nb"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2")
precYrs <- subset(prec, Coef=="beta" | Coef=="intYr")
years <- unique(allD$year)[2:14]+1900
precYrs$Year <- c(rep(rep(years, each=3000), each=4),
                   rep(rep(years, each=3000), each=4))

####
#### Now get subset defined by parameters; in a function
####
getRecCoefs <- function(doYear, mcDraw, doSpp){
  #First do coefficients with random year effects
  recNowYr <- subset(precYrs, Spp==doSpp & Year==doYear & Iter==mcDraw)
  iID <- which(recNowYr$Coef=="intYr")
  intercept <- recNowYr$value[iID]
  sID <- which(recNowYr$Coef=="beta")
  size <- recNowYr$value[sID]
  
  #Now do climate and competition fixed effects
  recNow <- subset(precAll, Spp==doSpp & Iter==mcDraw)
  cID <- which(recNow$Coef=="rain1"|recNow$Coef=="rain2"|recNow$Coef=="temp1"|recNow$Coef=="temp2")
  climEffs <- recNow$value[cID]
  dd <- recNow$value[which(recNow$Coef=="nb")]
  
  #Collate all parameters for output
  Rpars=list(intcpt=intercept, 
             slope=size,
             nb=dd,
             clim=climEffs)
  
  out(Rpars)
}
