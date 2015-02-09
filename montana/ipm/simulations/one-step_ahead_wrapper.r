
# this script loops through all observed year-to-year transitions,
# sets up initial conditions, and then calls the IPM to project the population

#--------------------------------------------------
# INPUTS
#--------------------------------------------------

#set working directory
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/montana/ipm/simulations",sep="")); # modify as needed 

doSpp<-"BOGR"
spp_list<-c("BOGR","HECO","PASM","POSE") # all Montana species
reps<-1   # number of times to simulate each quadrat x year transition
nMCMC<-3000 # max number of MCMC iterations to draw parameters from
outfile1<-paste(doSpp,"_sim_cover_1step_ahead.csv",sep="")
outfile2<-paste(doSpp,"_sim_density_1step_ahead.csv",sep="")

#--------------------------------------------------
# GET SIZE CLASSES
#--------------------------------------------------

sppCode<-which(spp_list==doSpp)
n_spp<-length(spp_list) # this is needed b/c all 4 spp parameters are imported at once

bigM=c(75,50,10,50)[sppCode]            #Set matrix dimension for each species
maxSize=c(2500,120,25,100)[sppCode]    

v=v.r=b.r=expv=Cr=WmatG=WmatS=NULL #a LIST for different species
h=r.L=r.U=Ctot=NULL

L=log(0.2)
U=log(maxSize)*1.1     
  
# b---boundary points. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
b = L+c(0:bigM)*(U-L)/bigM #log-transformed

# mid points
v = 0.5*(b[1:bigM]+b[2:(bigM+1)]) #log-transformed
 
#--------------------------------------------------
# GET OBSERVED DATA
#--------------------------------------------------

# Read in climate data 
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:5] <- scale(clim_data[2:5], center = TRUE, scale = TRUE) # standardize

# get calendar years
tmp <- read.csv("../../speciesData/quadAllCover.csv")
years <- unique(tmp$year)
years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year

# Read in quad inventory
quadYearList <- read.csv("../../speciesData/quad_inventory_2.csv")
quads<-names(quadYearList)

# Get quad-group info
quadGroups <- read.csv("../../speciesData/quad_info.csv")
names(quadGroups)<-c("quad","Group","grazing")

# Read in genet data
gen_dat<-read.csv(paste("../../speciesData/",doSpp,"/survD.csv",sep=""))

#--------------------------------------------------
# LOOP THROUGH QUADS AND YEARS
#--------------------------------------------------

output<-data.frame(quad="c",t0=1,t1=2,rep=1,cover.t0=1,cover.t1=1,stringsAsFactors=F)
counter<-0
for(iQ in 1:length(quads)){
  doGroup<-quadGroups$Group[which(quadGroups$quad==quads[iQ])]
  doGroup<-as.numeric(doGroup)  # from factor to scalar
  for(iYr in 1:length(years)){
    
    doYear<-years[iYr]
    weather<-clim_data[clim_data$year==(1900+doYear),2:5]
    
    # only work with complete transitions
    if(!is.na(quadYearList[iYr,iQ]) & !is.na(quadYearList[iYr+1,iQ])){
      
      tmpD<-subset(gen_dat,quad==quads[iQ] & year==years[iYr])
      nt.init<-table(cut(log(tmpD$area),breaks=b))
      cover.t0<-sum(tmpD$area)
       
      for(iRep in 1:reps){
          # call IPM script, passing doGroup, doYear, and nt.init
              cover.t1<-9999
          # store nt.new
          counter<-counter+1
          output[counter,]<-NA
          output$quad[counter]<-quads[iQ]
          output[counter,2:6]<-c(doYear,doYear+1,iRep,cover.t0,cover.t1)
      } # next iRep
       
    } # end if
     
    
  } # next iYr
} # next iQ 
