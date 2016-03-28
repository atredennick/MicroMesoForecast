##  Script to collate climate covariate mean and sds for scaling
##  Calculates scalars by species, vital rate, and left-out year

species <- c("BOGR", "HECO", "PASM", "POSE")
years <- c(1,33,34,35,36,37,38,39,40,41,42,43,44,45) #first year=1 for all years

full_df <- data.frame(species=NA, model=NA, covariate=NA, means=NA, sds=NA, yearout=NA)

####
####  PERCENT COVER COVARIATES
####
##  Read in data
#bring in data
allD <- read.csv("./speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("./weather/Climate.csv")
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


for(do_species in species){
  for(rm_year in 1:length(years)){
    growD <- subset(growD_all, Species==do_species & year!=years[rm_year])
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
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("qbm", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df[2:nrow(full_df),], "qbm_all_clim_scalars.RDS")



####
####  GROWTH PROCESS
####
full_df <- data.frame(species=NA, model=NA, covariate=NA, means=NA, sds=NA, yearout=NA)
sppList <- species
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
    sppD <- read.csv(paste("./speciesData/", doSpp, "/edited/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("./speciesData/", doSpp, "/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("./weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900

growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))
growD_all <- growD

for(do_species in species){
  for(rm_year in 1:length(years)){
    growD <- subset(growD_all, species==do_species & year!=years[rm_year])
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
    clim_covs <- growD[,clim_vars_all]
    # Get scalers for climate covariates from training data
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("growth", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df[2:nrow(full_df),], "growth_all_clim_scalars.RDS")



####
####  SURVIVAL PROCESS
####
full_df <- data.frame(species=NA, model=NA, covariate=NA, means=NA, sds=NA, yearout=NA)
sppList <- species
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

data_path <- "./speciesData/" #on local machine
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

climD <- read.csv("./weather/Climate.csv") #on local machine
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900

survD <- merge(survD,climD)
survD$Group=as.factor(substr(survD$quad,1,1))
survD_all <- survD

for(do_species in sppList){
  for(rm_year in 1:length(years)){
    survD <- subset(survD_all, species==do_species & year!=years[rm_year])
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
    clim_covs <- survD[,clim_vars_all]
    # Get scalers for climate covariates from training data
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("survival", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df[2:nrow(full_df),], "survival_all_clim_scalars.RDS")



####
####  RECRUITMENT PROCESS
####
climD <- read.csv("./weather/Climate.csv") #on local machine
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("./speciesData/", doSpp, "/edited/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("./speciesData/", doSpp, "/recArea.csv", sep=""))
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




###if using square transform, there is no need of the following code
####################################################################
###here you need to check: both parents1 and parents2 are equal to 0 at the same time
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

for(do_species in sppList){
  for(rm_year in 1:length(years)){
    recD <- subset(allD, species==paste("R.",do_species,sep="") & year!=years[rm_year])
    ##  Create and scale interaction covariates
    recD$ppt1TmeanSpr1 <- recD$ppt1*recD$TmeanSpr1
    recD$ppt2TmeanSpr2 <- recD$ppt2*recD$TmeanSpr2
    clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
    clim_covs <- recD[,clim_vars_all]
    clim_means <- colMeans(clim_covs)
    clim_sds <- apply(clim_covs, 2, FUN = sd)
    out_df <- data.frame(species = rep(do_species, length(clim_means)),
                         model = rep("recruitment", length(clim_means)),
                         covariate = names(clim_means),
                         means = as.numeric(clim_means),
                         sds = as.numeric(clim_sds),
                         yearout = years[rm_year])
    full_df <- rbind(full_df, out_df)
  }
}
saveRDS(full_df[2:nrow(full_df),], "recruitment_all_clim_scalars.RDS")