##  The script takes 'do_species' as a command line prompt. So, e.g.,
##    run as: "R CMD BATCH -1 survivalAllSpp_STAN.R" for species 1.

#clear everything, just to be safe 
# rm(list=ls(all=TRUE))

#load libraries
library(lme4)

sppList=sort(c("BOGR","HECO","PASM","POSE"))

####
#### Read in data by species and make one long data frame -------------
####
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area=NA,
                   survives=NA,
                   age=NA,
                   distEdgeMin=NA,
                   allEdge=NA,
                   species=NA)

data_path <- "../../../speciesData/" #on local machine
# data_path <- "speciesData/" #on HPC server

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste(data_path, doSpp, "/edited/survD.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste(data_path, doSpp, "/survD.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

survD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv") #on local machine
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900

survD <- merge(survD,climD)
survD$Group=as.factor(substr(survD$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRsurvCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOsurvCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMsurvCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEsurvCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)
colnames(crowd) <- c("W", "X", "species")

# Merge crowding and growth data
survD_all <- merge(survD, crowd, by=c("species", "X"))


##  Fit model to one species
survD <- subset(survD_all, species==do_species)
##  Create and scale interaction covariates
survD$ppt1TmeanSpr1 <- survD$ppt1*survD$TmeanSpr1
survD$ppt2TmeanSpr2 <- survD$ppt2*survD$TmeanSpr2
survD$sizepptLag <- survD$pptLag*log(survD$area)
survD$sizeppt1 <- survD$ppt1*log(survD$area)
survD$sizeppt2 <- survD$ppt2*log(survD$area)
survD$sizeTmeanSpr1 <- survD$TmeanSpr1*log(survD$area)
survD$sizeTmeanSpr2 <- survD$TmeanSpr2*log(survD$area)
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                   "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
survD[,clim_vars_all] <- scale(survD[,clim_vars_all], center = T, scale = T)
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
nyrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
survD$W2 <- survD$W*log(survD$area)
yid <- as.numeric(as.factor(survD$year))

surv_out1 <- glmer(survives~log(area)+W+W2+(1|Group)+(log(area)|year)+
               pptLag+ppt1+ppt2+TmeanSpr1+TmeanSpr2+ppt1TmeanSpr1+ppt2TmeanSpr2, 
             data=survD, family=binomial)

surv_out2 <- glmer(survives~log(area)+W+W2+(1|Group)+(log(area)|year)+
               pptLag+ppt1+ppt2+TmeanSpr1+TmeanSpr2+ppt1TmeanSpr1+ppt2TmeanSpr2+
               sizeppt1+sizeppt2+sizeTmeanSpr1+sizeTmeanSpr2, data=survD, family=binomial)
