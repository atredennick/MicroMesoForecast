# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Test script to pull in rect parameters for IPM #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#Several arguments are required
#' @param doYear  Specific climate year to pull year random effects (19xx)
#' @param mcDraw  A numeric scalar for the row(s) of MCMC to draw parameters from 
#' @param doSpp   A character scalar for the current speicies
#' @param group   A numeric scalar with group indicator

library(reshape2)
library(plyr)

####
#### Read in full MCMC output and format as data frame
####
MCMC <- readRDS("../vitalRateRegs/recruitment/recruitmentParamsMCMC.rds")
pRec2 <- melt(MCMC)
pRec2$Spp <- c(rep(spp_list, each=3000),
                rep(rep(spp_list, each=3000), times=6),
                rep(rep(spp_list, each=3000), times=13),
                rep(rep(spp_list, each=3000), times=6))
pRec2$Coef <- c(rep("dd", times=3000*4),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*13),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000),
                rep("theta", times=4*3000),
                rep("u", times=4*3000))
colnames(pRec2)[1] <- "Iter"
pRec <- pRec2[,c(1,3:5)]; rm(pRec2)
pRecAll <- subset(pRec, Coef=="gInt"|Coef=="dd"|Coef=="rain1"|Coef=="rain2"|Coef=="temp1"|Coef=="temp2"|Coef=="theta"|Coef=="u")
pRecYrs <- subset(pRec, Coef=="intYr")
pRecYrs$Year <- c(rep(rep(years, each=3000), each=4))

####
#### Now get subset defined by parameters; in a function
####
getRecCoefs <- function(doYear, mcDraw, group){
  #First do coefficients with random year effects
  recNowYr <- subset(pRecYrs, Year==doYear & Iter==mcDraw)
  iID <- which(recNowYr$Coef=="intYr")
  intercept <- recNowYr$value[iID]
  
  #Now do group, climate, and competition fixed effects
  recNow <- subset(pRecAll, Iter==mcDraw)
  cID <- which(recNow$Coef=="rain1"|recNow$Coef=="rain2"|recNow$Coef=="temp1"|recNow$Coef=="temp2")
  climEffs <- matrix(recNow$value[cID], 4, n_spp)
  dd <- recNow$value[which(recNow$Coef=="dd")]
  gID <- which(recNow$Coef=="gInt")
  ifelse(is.na(group)==TRUE,
         intG <- rep(0, n_spp),
         intG <- recNow$value[gID[group]])
  u <- recNow$value[which(recNow$Coef=="u")]
  theta <- recNow$value[which(recNow$Coef=="theta")]
  
  #Collate all parameters for output
  Rpars=list(intcpt=intercept, 
             grpInt=intG,
             dd=dd,
             clim=climEffs,
             u=u,
             theta=theta)
  return(Rpars)
}
