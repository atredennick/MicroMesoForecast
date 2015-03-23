##  This script runs validation models for the quad-based population model.

##  Statistical results are from leave-one-year-out fits, so here we attempt
##    to use the model results to predict the year left out of fitting. We
##    fit 13 different models, so this script loops through the 13 unfit years
##    to predict to cover value for that year. We also loop through 100 different
##    simulations for each quadrat in each year to capture parameter uncertainty 
##    by drawing parameter values from the MCMC chain for each simulation.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

##  Date began:     3.23.2015
##  Date completed: TBD



##  Clear the workspace
rm(list=ls(all=TRUE))

## Set do_year for validation
do_year <- 32

####
####  Load libraries -----------------------------------
####
library('reshape2')
library('plyr')
library('EnvStats')


####
####  Read in data -------------------------------------
####
##  Observations
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

##  Climate
climD <- read.csv("../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)


####
####  Read in statistical model parameters ------------
####
in_file <- paste("../vitalRateRegressions/truncNormModel/validation/growthParamsMCMC_", do_year, ".rds", sep="")
pGrow <- readRDS(in_file)
pGrow2 <- melt(pGrow)
pGrow2$Spp <- c(rep(rep(sppList, each=3000), times=12),
                rep(rep(sppList, each=3000), times=1),
                rep(rep(sppList, each=3000), times=6),
                rep(rep(sppList, each=3000), times=12),
                rep(rep(sppList, each=3000), times=7))
pGrow2$Coef <- c(rep("beta", times=3000*4*12),
                 rep("betaSpp", times=3000*4),
                 rep("gInt", times=6*4*3000),
                 rep("intYr", times=4*3000*12),
                 rep("intercept", times=3000*4),
                 rep("rain1", times=4*3000),
                 rep("rain2", times=4*3000),
                 rep("rainLag", times=4*3000),
                 rep("tau", times=4*3000),
                 rep("temp1", times=4*3000),
                 rep("temp2", times=4*3000))
colnames(pGrow2)[1] <- "Iter"
pGrow <- pGrow2[,c(1,3:5)]; rm(pGrow2)
pGrowAll <- subset(pGrow, Coef=="gInt"|Coef=="rain1"|
                          Coef=="rain2"|Coef=="rainLag"|
                          Coef=="temp1"|Coef=="temp2"|
                          Coef=="tau"|Coef=="betaSpp"|Coef=="intercept")

