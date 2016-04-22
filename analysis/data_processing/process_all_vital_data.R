##  R Script to Process and Transform Raw Data into Dataframes Used
##  for Statistical Modeling and Model Simulations.
##
##  Author: Andrew Tredennick
##  Last update: 4-21-2016

# Clear the workspace
rm(list=ls())


####
####  LIBRARIES
####
library("plyr")
library("reshape2")



####
####  PRELIMINARIES
####
species <- c("BOGR", "HECO", "PASM", "POSE")
path2data <- "./speciesData/"
outpath <- "../processed_data/"
surv_outfile <- "surv_with_weather.RDS"
grow_outfile <- "grow_with_weather.RDS"
rec_outfile <- "rec_with_weather.RDS"
cover_outfile <- "cover_with_weather.RDS"
weather <- read.csv("../weather/Climate.csv")
my_covariates <- c("year","pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
weather <- weather[,my_covariates]



####
####  SOURCE FILES
####
source("dataExtractQuad.R") # Adds in zeros for all measured quads; rms bad bogr quadyrs
source("format_quad_data_fxn.R") # Function for formatting quadrat data



####
####  QUADRAT PROPORTIONAL COVER DATA
####
quad_dat <- read.csv(paste0(path2data,"quadAllCover.csv"))
out_quad_dat <- format_quad_data(quad_dat, weather)

# Add weather interactions
out_quad_dat$ppt1TmeanSpr1 <- with(out_quad_dat, ppt1*TmeanSpr1)
out_quad_dat$ppt2TmeanSpr2 <- with(out_quad_dat, ppt2*TmeanSpr2)

# Save
saveRDS(out_quad_dat, paste0(outpath,cover_outfile))



####
####  SURVIVAL DATA
####
surv_data <- list() # empty object to collate survival data frames
for(do_species in species){
  if(do_species=="BOGR") { getfile <- paste0(path2data,do_species,"/edited/survD.csv") }
  if(do_species!="BOGR") { getfile <- paste0(path2data,do_species,"/survD.csv") }
  tmpsurv <- read.csv(getfile)
  tmpsurv$species <- do_species
  surv_data <- rbind(surv_data, tmpsurv)
}

# Get rid of the 1900 for year in a weather df copy; to match survival df year
weather_copy <- weather
weather_copy$year <- weather_copy$year-1900

# Merge in weather data
surv_data <- merge(surv_data,weather_copy)

# Extract group information
surv_data$Group <- as.factor(substr(surv_data$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRsurvCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOsurvCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMsurvCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEsurvCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowdS <- rbind(c1,c2,c3,c4)
colnames(crowdS) <- c("W", "X", "species")

# Merge crowding and growth data
surv_data <- merge(surv_data, crowdS, by=c("species", "X"))
# Add weather interactions
surv_data$ppt1TmeanSpr1 <- with(surv_data, ppt1*TmeanSpr1)
surv_data$ppt2TmeanSpr2 <- with(surv_data, ppt2*TmeanSpr2)

# Save the survival data frame
saveRDS(surv_data, paste0(outpath,surv_outfile))



####
####  GROWTH DATA
####
growth_data <- list() # empty object to collate growth data frames
for(do_species in species){
  if(do_species=="BOGR") { getfile <- paste0(path2data,do_species,"/edited/growDnoNA.csv") }
  if(do_species!="BOGR") { getfile <- paste0(path2data,do_species,"/growDnoNA.csv") }
  tmpgrowth <- read.csv(getfile)
  tmpgrowth$species <- do_species
  growth_data <- rbind(growth_data, tmpgrowth)
}

# Get rid of the 1900 for year in a weather df copy; to match growth df year
weather_copy <- weather
weather_copy$year <- weather_copy$year-1900

# Merge in weather data
growth_data <- merge(growth_data,weather_copy)

# Extract group ID information
growth_data$Group <- as.factor(substr(growth_data$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRgrowthCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOgrowthCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMgrowthCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEgrowthCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowdG <- rbind(c1,c2,c3,c4)
colnames(crowdG) <- c("W", "X", "species")

# Merge crowding and growth data
growth_data <- merge(growth_data, crowdG, by=c("species", "X"))

# Add weather interactions
growth_data$ppt1TmeanSpr1 <- with(growth_data, ppt1*TmeanSpr1)
growth_data$ppt2TmeanSpr2 <- with(growth_data, ppt2*TmeanSpr2)

# Save the growth data frame
saveRDS(growth_data, paste0(outpath,grow_outfile))



####
####  RECRUITMENT DATA
####
for(do_species in species){
  if(do_species=="BOGR") { getfile <- paste0(path2data,do_species,"/edited/recArea.csv") }
  if(do_species!="BOGR") { getfile <- paste0(path2data,do_species,"/recArea.csv") }
  tmprec <- read.csv(getfile)
  tmprec$species <- do_species
  tmprec$Group <- as.factor(substr(tmprec$quad,1,1)) # extract group information
  tmprec <- tmprec[,c("quad","year","NRquad","totParea","Group")]
  names(tmprec)[3] <- paste0("R.",do_species)
  names(tmprec)[4] <- paste0("cov.",do_species)
  
  sppid <- which(species==do_species)
  if(sppid==1){
    rec_data <- tmprec
  }else{
    rec_data <- merge(rec_data,tmprec,all=T)
  }
}
rec_data[is.na(rec_data)]=0  # replace missing values 

D <- rec_data # rename for easy merging with Chengjin's code...

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
rec_final <- data.frame(species=tmpY$Var2,
                   year=rep(year,4),
                   group=rep(Group,4),
                   recruits=tmpY$value,
                   parents1=tmpP1$value,
                   parents2=tmpP2$value)


# Get rid of the 1900 for year in a weather df copy; to match growth df year
weather_copy <- weather
weather_copy$year <- weather_copy$year-1900
rec_final <- merge(rec_final,weather_copy)

# Add weather interactions
rec_final$ppt1TmeanSpr1 <- with(rec_final, ppt1*TmeanSpr1)
rec_final$ppt2TmeanSpr2 <- with(rec_final, ppt2*TmeanSpr2)

# Save the growth data frame
saveRDS(rec_final, paste0(outpath,rec_outfile))

