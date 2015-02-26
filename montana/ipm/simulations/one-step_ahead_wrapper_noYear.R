
# this script loops through all observed year-to-year transitions,
# sets up initial conditions, and then calls the IPM to project the population

#--------------------------------------------------
# INPUTS
#--------------------------------------------------

#set working directory
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/montana/ipm/simulations",sep="")); # modify as needed 

# doSpp<-"BOGR"
spp_list<-c("BOGR","HECO","PASM","POSE") # all Montana species
reps<-100  # number of times to simulate each quadrat x year transition
nMCMC<-3000 # max number of MCMC iterations to draw parameters from

#outfile2<-paste(doSpp,"_sim_density_1step_ahead.csv",sep="") # not implemented yet

#--------------------------------------------------
#  Set up simulation, Matrix size and initial vectors
#--------------------------------------------------
for(spp in 1:length(spp_list)){
  doSpp <- spp_list[spp]
  outfile1<-paste("./one_step_results/",doSpp,"_sim_cover_1step_ahead_noYear.csv",sep="")
  sppCode<-which(spp_list==doSpp)
  n_spp<-length(spp_list) # this is needed b/c all 4 spp parameters are imported at once
  
  Atotal<-10000 #area of quadrat in cm^2
  
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
  
  # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
  h = v[2]-v[1] #log-transformed 
  
  # variables for Wr approximation---radius
  b.r=sqrt(exp(b)/pi) 
  v.r=sqrt(exp(v)/pi)
  expv=exp(v) #for size
  r.L=sqrt(exp(L)/pi) #the lower size limit; radius
  r.U=sqrt(exp(U)/pi) #the upper size limit; radius
  WmatG=rep(NA,length(v.r))  # storage of size-specific conspecific W values for each species
  WmatS=rep(NA,length(v.r))  # storage of size-specific conspecific W values for each species
  
  tmp=range(v.r)#this is the 'real" range, not log-transformed...
  size.range=seq(tmp[1],tmp[2],length=bigM) # range across all possible sizes; 'real' radius for each size stage
  
  # read in IPM functions
  source("one_step_ahead_functions.r")
  
  #--------------------------------------
  # VITAL RATE FUNCTIONS AND PARAMETERS
  #--------------------------------------
  
  source("vital_rate_ipm_functions.R")
  
  # get calendar years
  raw_data <- read.csv("../../speciesData/quadAllCover.csv")
  years <- unique(raw_data$year)
  years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year
  
  source("../vitalRateRegs/survival/import2ipm.R")
  source("../vitalRateRegs/growth/import2ipm.R")
  source("../vitalRateRegs/recruitment/import2ipm.R")
  
  # get recruit size parameters
  rec_size_mean <- numeric(n_spp)
  rec_size_var <- numeric(n_spp)
  for(i in 1:n_spp){
    infile=paste("../../speciesData/",spp_list[i],"/recSize.csv",sep="")
    recSize=read.csv(infile)
    rec_size_mean[i]=mean(log(recSize$area))
    rec_size_var[i]=var(log(recSize$area))
  }
  
  # get alphas values (needed to calculate neighborhood crowding)
  alpha_grow <- read.csv("../../alpha_list_growth.csv")
  alpha_surv <- read.csv("../../alpha_list_survival.csv")
  # pull out alphas only for the doSpp
  alphaG <- subset(alpha_grow, Site=="Montana")$Alpha[sppCode]
  alphaS <- subset(alpha_surv, Site=="Montana")$Alpha[sppCode]
  
  #--------------------------------------------------
  # GET OBSERVED DATA
  #--------------------------------------------------
  
  # Read in climate data 
  clim_data <- read.csv("../../weather/Climate.csv")
  clim_data <- clim_data[,c("year","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
  clim_data[2:5] <- scale(clim_data[2:5], center = TRUE, scale = TRUE) # standardize
  
  # get calendar years
  coverDat <- read.csv("../../speciesData/quadAllCover.csv")
  years <- unique(coverDat$year)
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
      
      climYear<-years[iYr]
      weather<-clim_data[clim_data$year==(1900+climYear),2:5]
      doYear<-NA
      
      # only work with complete transitions
      if(!is.na(quadYearList[iYr,iQ]) & !is.na(quadYearList[iYr+1,iQ])){
        
        tmpD<-subset(gen_dat,quad==quads[iQ] & year==years[iYr])
        nt.init<-table(cut(log(tmpD$area),breaks=b))
        cover.t0<-sum(tmpD$area)/Atotal
        
        if(cover.t0>0){
          
          # create Cr function for calculating neighborhood crowding
          # (I don't understand why this doesn't work when it is inside projectIPM() )
          #Ctot=h*sum(expv*nt.init) #h*sum(expv*nt.init) #total cover--multiply by h or not?
          #Cr=splinefun(b.r,h*c(0,cumsum(expv*nt.init)),method="natural") #Cr is a function             
          Ctot=sum(expv*nt.init) #h*sum(expv*nt.init) #total cover--multiply by h or not?
          Cr=splinefun(b.r,c(0,cumsum(expv*nt.init)),method="natural") #Cr is a function             
          
          
          for(iRep in 1:reps){
            
            # parameter draw
            mcDraw<-sample(1:nMCMC,1) 
            
            # call IPM script
            nt.new<-projectIPM(nt=nt.init,doYear,doGroup,mcDraw,weather,sppCode)
            cover.t1<-sum(nt.new*exp(v))/Atotal
            
            # store nt.new
            counter<-counter+1
            output[counter,]<-NA
            output$quad[counter]<-quads[iQ]
            output[counter,2:6]<-c(doYear,doYear+1,iRep,cover.t0,cover.t1)
            
          } # next iRep
          
        } # end cover.t0>0 if
        
      } # end is.na(year) if
      
      print(paste(quads[iQ],years[iYr],sep=" "))
      flush.console()
      
    } # next iYr
    
  } # next iQ 
  
  # merge with observed data
  coverDat<-subset(coverDat,Species==doSpp)
  coverDat<-coverDat[,c("quad","year","propCover")]
  names(coverDat)<-c("quad","t1","obs.cover.t1")
  output<-merge(output,coverDat)
  
  plot(output$obs.cover.t1,output$cover.t1)
  abline(0,1)
  
  # write output
  write.csv(output,outfile1,row.names=F)
}




