##  R Script to collate summaries of data used (e.g., number of quadrats per
##  species, number of genets per species per year, etc.)

# Clear the workspace
rm(list=ls())

library("plyr")

####
####  PRELIMINARIES
####
species <- c("BOGR", "HECO", "PASM", "POSE")
quadfile <- "./speciesData/quadAllCover.csv"
path2gendat <- "./speciesData/"
weatherfile <- "../weather/Climate.csv"
source("format_quad_data_fxn.R")
"%w/o%" <- function(x, y) x[!x %in% y] # x without y


####
####  LOOK AT QUADRAT LEVEL DATA
####
quad_dat <- read.csv(quadfile) # quadrat data
weather_dat <- read.csv(weatherfile) # weather data
quad4mod <- format_quad_data(quad_dat,weather_dat) # formats data frame for modeling

# Summarise data by species
quad_Ns <- ddply(quad4mod, .(species), summarise,
                 N = length(propCover.t0),
                 Nq = length(unique(quad)))
# Store quadrats by species
quad_quads <- ddply(quad4mod, .(species), summarise,
                    quads = unique(quad))



####
####  LOOK AT GENET SURVIVAL DATA
####
gen_survdat <- list() # empty object that transforms to df on the fly
for(spp in species){
  if(spp=="BOGR") { gen_dat <- read.csv(paste0(path2gendat,spp,"/edited/survD.csv")) }
  if(spp!="BOGR") { gen_dat <- read.csv(paste0(path2gendat,spp,"/survD.csv")) }
  gen_dat <- gen_dat[ ,which(colnames(gen_dat)!="X")]
  gen_dat$species <- spp
  gen_survdat <- rbind(gen_survdat, gen_dat)
}
surv_Ns <- ddply(gen_survdat, .(species), summarise,
                 N = length(area),
                 Nq = length(unique(quad)))
# Store quadrats by species
surv_quads <- ddply(gen_survdat, .(species), summarise,
                    quads = unique(quad))

extra_quads <- as.character(surv_quads[which(surv_quads$species=="BOGR"),"quads"]) %w/o% quad_quads[which(quad_quads$Species=="BOGR"),"quads"] 
extra_survdat <- subset(gen_survdat, species=="BOGR" & quad%in%extra_quads)
gen_survdat[which(gen_survdat$trackID==-Inf),]

tmp <- ddply(gen_survdat, .(species, year, quad), summarise,
                         cover = sum(area/10000))
gen_based_cover <- ddply(tmp, .(species), summarise,
                         avg_cover = mean(cover)*100)
quad_based_cover <- ddply(quad4mod, .(species), summarise,
                          avg_cover = mean(propCover.t0)*100)

tmp$year <- tmp$year+1900
gen_quad <- merge(tmp,quad4mod, all.x=TRUE)

# remove NA rows where extirpition or colonization happend in the quad DF
Nrms_due2_extcol <- length(which(complete.cases(gen_quad)==FALSE))
gen_quad <- gen_quad[complete.cases(gen_quad),]
gen_quad_comparison <- ddply(gen_quad, .(species), summarise,
                             avg_gencov = mean(cover),
                             avg_quadcov = mean(propCover.t0))



####
####  MAKE SURE DISCREPENCY BETWEEN GENET DATA AND QUAD DATA COMES FROM
####  ELIMINATING ZEROS IN THE QUADRAT DATA FRAME
####
# remake quadrat data from, not using the function
quad_df <- quad_dat <- read.csv(quadfile) # quadrat data
climate_df <- weather_dat <- read.csv(weatherfile) # weather data
backD <- data.frame(climYear=NA,quad=NA,year=NA,
                    totCover=NA,Species=NA,propCover=NA,group=NA,
                    lag_propCover=NA,
                    pptLag=NA,ppt1=NA,TmeanSpr1=NA,ppt2=NA,TmeanSpr2=NA,
                    TmeanSum1=NA,TmeanSum2=NA)

sppList <- unique(quad_df[,"Species"])

# loop through species and remake data frame
for(ispp in 1:length(sppList)){
  doSpp <- sppList[ispp]
  sppD <- subset(quad_df, Species==doSpp)
  sppD$group <- substring(sppD$quad, 1, 1)
  
  # create lag cover variable
  tmp <- sppD[,c("quad","year","propCover")]
  tmp$year <- tmp$year+1
  names(tmp)[3] <- "lag_propCover"
  sppD <- merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear <- sppD$year+1900-1 # climate to be associated with lag cover (t0) 
  sppD <- merge(sppD,climate_df,by.x="climYear",by.y="year",keep.all=T)
  
  # Subset growth observations (no colonization or extirpation)
  # growD <- subset(sppD,lag_propCover>0 & propCover>0)
  growD <- sppD
  backD <- rbind(backD, growD)
}#end species loop

backD <- backD[complete.cases(backD),] # rms NA row
quadD <- data.frame(year=backD$climYear, 
                   species=backD$Species,
                   quad=backD$quad,
                   group=backD$group,
                   propCover.t1=backD$propCover,
                   propCover.t0=backD$lag_propCover)

tmp <- ddply(gen_survdat, .(species, year, quad), summarise,
             cover = sum(area/10000))
tmp$year <- tmp$year+1900
gen_quad2 <- merge(tmp,quadD, all.x=TRUE, all.y=TRUE)

##  THESE SHOULD BE IN THE BAD BOGR QUAD-YEARS 
##  (where the year in gen_quad is year-1 in the bad_bogr spreadsheet)
##  So, here the deal is that the survival data frame will have values for a quad-year one year
##  before a bad bogr quad-year (e.g., removed quad-year) but the quadrat data will not because
##  it is like the growth data frame and does not include rows where t1 or t0 are NA.
gen_quad2[which(is.na(gen_quad2$propCover.t0)==TRUE),c("species","year","quad")]


## Just look at the raw data frames
genet_based_cover <- tmp
quad_based_cover <- read.csv("speciesData/quadAllCover.csv")
names(quad_based_cover) <- tolower(names(quad_based_cover))
quad_based_cover$year <- quad_based_cover$year+1900
all <- merge(genet_based_cover,quad_based_cover,all.x=T,all.y=T)
all_noNA_noZeros <- all[which(is.na(all$cover)==FALSE),] #rms NAs and zeros b/c quad df has zeros where genet df has NA once merged
avg_cov <- colMeans(all_noNA_noZeros[,c("cover","propcover")])
if(avg_cov[1]!=avg_cov[2]){ stop("YOU HAVE A PROBLEM HERE") }