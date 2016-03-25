##
##  R script for survival logistic model as run on Utah State University's
##  High Performance Computing system.
##
##  Script takes command line argument 'do_grid' to run all levels of
##  regularization and each cross-validation set in parallel. Script launches
##  n*k models: 'n' prior standard deviations and 'k' validation sets.
##
##  Saved output:
##    1. Summed log pointwise predictive density based on 
##       out-of-sample predictions
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 1-20-2016
##

# Change species four letter code here...
do_species <- "BOGR"

####
####  Set SD Prior and CV Set from Command Line Arguments
####
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_grid <- as.numeric(myargument)


####
####  Load Libraries and Subroutines
####
library(rstan)
library(plyr)
library(reshape2)
library(ggmcmc)
library(matrixStats)
source("waic_fxns.R")



####
####  Load Climate Covariates and Recruitment Data
####
climD <- read.csv("Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900


sppList=sort(c("BOGR","HECO","PASM","POSE"))

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("./", doSpp, "/edited/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("./", doSpp, "/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }
  
  sppD$Group <- as.factor(substr(sppD$quad,1,1)) #add by Chengjin
  sppD <- sppD[,c("quad","year","NRquad","totParea","Group")]
  names(sppD)[3] <- paste("R.",sppList[spp],sep="")
  names(sppD)[4] <- paste("cov.",sppList[spp],sep="")
  if(spp==1){
    D <- sppD
  }else{
    D <- merge(D,sppD,all=T)
  }
}
D[is.na(D)]=0  # replace missing values 

# calculate mean cover by group and year
tmpD <- D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD <- aggregate(tmpD[,4:NCOL(tmpD)],
                  by=list("year"=tmpD$year,"Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)] <- paste("Gcov.",sppList,sep="")
D <- merge(D,tmpD,all.x=T)

parents1 <- as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100 ##convert from absolute cover to [1,100] range
parents2 <- as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100

##for species 1
tmp1L <- which(parents1[,1]==0) ##lcoal
tmp1G <- which(parents2[,1]==0) ##Group
tmp1 <- intersect(tmp1L,tmp1G)
##for species 2
tmp2L <- which(parents1[,2]==0)
tmp2G <- which(parents2[,2]==0)
tmp2 <- intersect(tmp2L,tmp2G)
##for species 3
tmp3L <- which(parents1[,3]==0)
tmp3G <- which(parents2[,3]==0)
tmp3 <- intersect(tmp3L,tmp3G)
##for species 4
tmp4L <- which(parents1[,4]==0)
tmp4G <- which(parents2[,4]==0)
tmp4 <- intersect(tmp4L,tmp4G)

tmp <- unique(c(tmp1,tmp2,tmp3,tmp4))

if(length(tmp)>0){
  parents1 <- parents1[-tmp,] ##remove them
  parents2 <-parents2[-tmp,] ##remove them
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])[-tmp,] ##remove them  
  year <- D$year[-tmp] ##remove them
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]-length(tmp) ##reduce
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
} else {
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])
  year <- D$year
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
}

tmpY <- melt(y)
tmpP1 <- melt(parents1)
tmpP2 <- melt(parents2)
allD <- data.frame(species=tmpY$Var2,
                   year=rep(year,4),
                   group=rep(Group,4),
                   recruits=tmpY$value,
                   parents1=tmpP1$value,
                   parents2=tmpP2$value)

##  Add in climate data
allD <- merge(allD,climD)



####
####  Compile Stan Model
####
holdyear <- unique(allD$year)[1]
trainD <- subset(allD, species==paste0("R.",do_species) & year!=holdyear)
##  Create and scale interaction covariates
trainD$ppt1TmeanSpr1 <- trainD$ppt1*trainD$TmeanSpr1
trainD$ppt2TmeanSpr2 <- trainD$ppt2*trainD$TmeanSpr2
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- trainD[,clim_vars_all]

groups <- as.numeric(trainD$group)
G <- length(unique(trainD$group))
Yrs <- length(unique(trainD$year))
yid <- as.numeric(as.factor(trainD$year))

# Get scalers for climate covariates from training data
clim_means <- colMeans(clim_covs)
clim_sds <- apply(clim_covs, 2, FUN = sd)
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)

holdD <- subset(allD, species==paste0("R.",do_species) & year==holdyear)
holdD$ppt1TmeanSpr1 <- holdD$ppt1*holdD$TmeanSpr1
holdD$ppt2TmeanSpr2 <- holdD$ppt2*holdD$TmeanSpr2
clim_covs_oos <- holdD[,clim_vars_all]
for(j in 1:ncol(clim_covs_oos)){
  clim_covs_oos[,j] <- (clim_covs_oos[,j] - clim_means[j])/clim_sds[j]
}
groups_out <- as.numeric(holdD$group)
Gout <- length(unique(holdD$group))
npreds <- nrow(holdD)
y_holdout <- holdD$recruits
parents1_out <- holdD$parents1
parents2_out <- holdD$parents2

datalist <- list(N=nrow(trainD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=trainD$recruits, C=clim_covs, 
                 parents1=trainD$parents1, parents2=trainD$parents2,
                 G=G, gid=groups, tau=1,
                 npreds=npreds, y_holdout=y_holdout, C_out=clim_covs_oos,
                 parents1_out=parents1_out, parents2_out=parents2_out,
                 gid_out=groups_out, Gout=Gout)

pars=c("log_lik")
mcmc_oos <- stan(file = "recruitment_oos_cv.stan", data=datalist,
                 pars=pars, chains=0)



####
####  Set up CV x Regularization grid
####
recD <- subset(allD, species==paste0("R.",do_species))
n.beta <- 24
sd_vec <- seq(0.1,1.5,length.out = n.beta)
yrs.vec <- unique(allD$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)



####
####  Source Recruitment Model Function and Fit Models
####
source("recruitment_fxns.R")
out_lpd <- cv.fcn(do_grid)
saveRDS(out_lpd, paste0(do_species,"_oos_cv_dogrid_",do_grid,".RDS"))

