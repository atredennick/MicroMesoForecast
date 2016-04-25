##  Function to transform "raw" quadrat data into population growth
##  dataframe for model fitting.

format_quad_data <- function(quad_df, climate_df){
  backD <- list()
  
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
    sppD <- merge(sppD,climate_df,by.x="climYear",by.y="year")
    
    # Subset growth observations (no colonization or extirpation)
    growD <- subset(sppD,lag_propCover>0 & propCover>0)
    backD <- rbind(backD, growD)
  }#end species loop
  
  backD <- backD[complete.cases(backD),] # rms NA row
  outD <- data.frame(year=backD$climYear, 
                     species=backD$Species,
                     quad=backD$quad,
                     group=backD$group,
                     propCover.t1=backD$propCover,
                     propCover.t0=backD$lag_propCover,
                     pptLag=backD$pptLag,
                     ppt1=backD$ppt1,
                     ppt2=backD$ppt2,
                     TmeanSpr1=backD$TmeanSpr1,
                     TmeanSpr2=backD$TmeanSpr2)
  return(outD)
}