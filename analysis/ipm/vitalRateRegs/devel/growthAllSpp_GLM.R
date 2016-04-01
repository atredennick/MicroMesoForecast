#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
# rm(list=ls(all=TRUE))

#load libraries
library(lme4)
sppList=sort(c("BOGR","HECO","PASM","POSE"))
do_species <- "POSE"
####
#### Read in data by species and make one long data frame -------------
####
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area.t1=NA,
                   area.t0=NA,
                   age=NA,
                   allEdge=NA,
                   distEdgeMin=NA,
                   species=NA)

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900

growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRgrowthCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOgrowthCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMgrowthCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEgrowthCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)

## Merge crowding and growth data
growD_all <- merge(growD, crowd, by=c("species", "X"))

## Compile model outside of loop
growD <- subset(growD_all, species==do_species)
##  Create and scale interaction covariates
growD$ppt1TmeanSpr1 <- growD$ppt1*growD$TmeanSpr1
growD$ppt2TmeanSpr2 <- growD$ppt2*growD$TmeanSpr2
growD$sizepptLag <- growD$pptLag*log(growD$area.t0)
growD$sizeppt1 <- growD$ppt1*log(growD$area.t0)
growD$sizeppt2 <- growD$ppt2*log(growD$area.t0)
growD$sizeTmeanSpr1 <- growD$TmeanSpr1*log(growD$area.t0)
growD$sizeTmeanSpr2 <- growD$TmeanSpr2*log(growD$area.t0)
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                   "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
growD[,clim_vars_all] <- scale(growD[,clim_vars_all], center = T, scale = T)
# clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(growD$Group)
G <- length(unique(growD$Group))
nyrs <- length(unique(growD$year))
W <- cbind(growD$W, growD$W*log(growD$area.t0))
yid <- as.numeric(as.factor(growD$year))

growth_out1 <- lmer(log(area.t1)~log(area.t0)+W+(1|Group)+(log(area.t0)|year)+
              pptLag+ppt1+ppt2+TmeanSpr1+TmeanSpr2+ppt1TmeanSpr1+ppt2TmeanSpr2, 
             data=growD)

growth_out2 <- lmer(log(area.t1)~log(area.t0)+W+(1|Group)+(log(area.t0)|year)+
               pptLag+ppt1+ppt2+TmeanSpr1+TmeanSpr2+ppt1TmeanSpr1+ppt2TmeanSpr2+
               sizepptLag+sizeppt1+sizeppt2+sizeTmeanSpr1+sizeTmeanSpr2, data=growD)

AIC(growth_out1,growth_out2)

