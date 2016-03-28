## Script to calculate species-specific climate covariate scalars

rm(list=ls())

##  Read in data
#bring in data
allD <- read.csv("../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")

backD <- data.frame(climYear=NA,
                    quad = NA,
                    year= NA,
                    totCover= NA,
                    Species= NA,
                    propCover= NA,
                    lag.cover= NA,
                    pptLag= NA,
                    ppt1= NA,
                    TmeanSpr1= NA,
                    ppt2= NA,
                    TmeanSpr2= NA,
                    TmeanSum1= NA,
                    TmeanSum2= NA,
                    yearID= NA,
                    group = NA,
                    percCover = NA,
                    percLagCover = NA)

#loop through species and remake data frame
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover variable
  tmp=sppD[,c("quad","year","totCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000
  backD <- rbind(backD, growD)
}#end species loop
growD_all <- backD[2:nrow(backD),]


for(do_species in sppList){
  growD <- subset(growD_all, Species=="BOGR")
  ##  Create and scale interaction covariates
  growD$ppt1TmeanSpr1 <- growD$ppt1*growD$TmeanSpr1
  growD$ppt2TmeanSpr2 <- growD$ppt2*growD$TmeanSpr2
  growD$sizepptLag <- growD$pptLag*log(growD$percLagCover)
  growD$sizeppt1 <- growD$ppt1*log(growD$percLagCover)
  growD$sizeppt2 <- growD$ppt2*log(growD$percLagCover)
  growD$sizeTmeanSpr1 <- growD$TmeanSpr1*log(growD$percLagCover)
  growD$sizeTmeanSpr2 <- growD$TmeanSpr2*log(growD$percLagCover)
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2", "sizepptLag",
                     "sizeppt1", "sizeppt2", "sizeTmeanSpr1", "sizeTmeanSpr2")
  clim_covs <- growD[,clim_vars_all]
  # Get scalers for climate covariates from training data
  clim_means <- colMeans(clim_covs)
  clim_sds <- apply(clim_covs, 2, FUN = sd)
  out_df <- data.frame(covariate = names(clim_means),
                       means = clim_means,
                       sds = clim_sds)
  saveRDS(out_df, paste0(do_species, "climscalars.RDS"))
}
