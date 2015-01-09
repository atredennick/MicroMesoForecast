# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Test script to pull in growth parameters for IPM #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#Several arguments are required
#' @param doYear  Specific climate year to pull year random effects
#' @param mcDraw  A numeric scalar or vector for the row(s) of MCMC to draw parameters from 
#' @param doSpp   A character scalar for the current speicies
#' @param group   A numeric scalar with group indicator

library(reshape2)
library(plyr)

####
#### Read in full MCMC output and format as data frame
####
MCMC <- readRDS("growthParamsMCMC.rds")
pGrow2 <- melt(MCMC)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=13),
                rep(rep(sppList, each=3000), times=5))
pGrow2$Coef <- c(rep("beta", times=3000*4*13),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*13),
                 rep("nb", times=4*3000),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(pGrow2)[1] <- "Iter"
pGrow <- pGrow2[,c(1,3:5)]; rm(pGrow2)
pGrowAll <- subset(pGrow, Coef=="gInt"|Coef=="nb"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2")
pGrowYrs <- subset(pGrow, Coef=="beta" | Coef=="intYr")
pGrowYrs$Year <- c(rep(rep(years, each=3000), each=4),
                   rep(rep(years, each=3000), each=4))

####
#### Now get subset defined by parameters; in a function
####
getGrowCoefs <- function(doYear, mcDraw, doSpp){
  #First do coefficients with random year effects
  growNowYr <- subset(pGrowYrs, Spp==doSpp & Year==doYear & Iter==mcDraw)
  iID <- which(growNowYr$Coef=="intYr")
  intercept <- growNowYr$value[iID]
  sID <- which(growNowYr$Coef=="beta")
  size <- growNowYr$value[sID]
  
  #Now do climate and competition fixed effects
  growNow <- subset(pGrowAll, Spp==doSpp & Iter==mcDraw)
  cID <- which(growNow$Coef=="rain1"|growNow$Coef=="rain2"|growNow$Coef=="temp1"|growNow$Coef=="temp2")
  climEffs <- growNow$value[cID]
  dd <- growNow$value[which(growNow$Coef=="nb")]
  gID <- which(growNow$Coef=="gInt")
  intG <- growNow$value[gID[group]]
  
  #Get variance parameters
  varPars <- readRDS("varianceParams.rds")
  
  #Collate all parameters for output
  Gpars=list(intcpt=intercept, 
             intG=intG,
             slope=size,
             nb=dd,
             clim=climEffs,
             sigma2.a=subset(varPars, species==doSpp)$a,
             sigma2.b=subset(varPars, species==doSpp)$b)
  
  out(Gpars)
}
