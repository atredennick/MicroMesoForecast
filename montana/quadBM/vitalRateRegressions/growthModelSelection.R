# IBM version: each "individual" is a quadrat, 
# we track cover in each quadrat over time

#This script is for model selection based on DIC, regardless of 95% CIs

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

####
#### POSSIBLE FIXED-EFFECTS -----------------------
####
fixed.forms <- list(form1 = "log(totCover) ~ log(lag.cover)+ppt1*TmeanSpr1+ppt2*TmeanSpr2",
                    form2 = "log(totCover) ~ log(lag.cover)+ppt1*TmeanSpr1+ppt2+TmeanSpr2",
                    form3 = "log(totCover) ~ log(lag.cover)+ppt1+TmeanSpr1+ppt2*TmeanSpr2",
                    form4 = "log(totCover) ~ log(lag.cover)+ppt1*TmeanSpr1+ppt2",
                    form5 = "log(totCover) ~ log(lag.cover)+ppt1*TmeanSpr1+TmeanSpr2",
                    form6 = "log(totCover) ~ log(lag.cover)+ppt1+ppt2*TmeanSpr2",
                    form7 = "log(totCover) ~ log(lag.cover)+TmeanSpr1+ppt2*TmeanSpr2",
                    form8 = "log(totCover) ~ log(lag.cover)+ppt1+TmeanSpr1+ppt2+TmeanSpr2",
                    form9 = "log(totCover) ~ log(lag.cover)+ppt1+TmeanSpr1+ppt2",
                    form10 = "log(totCover) ~ log(lag.cover)+ppt1+TmeanSpr1+TmeanSpr2",
                    form11 = "log(totCover) ~ log(lag.cover)+ppt1+ppt2+TmeanSpr2",
                    form12 = "log(totCover) ~ log(lag.cover)+TmeanSpr1+ppt2+TmeanSpr2",
                    form13 = "log(totCover) ~ log(lag.cover)+ppt1+ppt2",
                    form14 = "log(totCover) ~ log(lag.cover)+ppt1+TmeanSpr2",
                    form15 = "log(totCover) ~ log(lag.cover)+TmeanSpr1+ppt2",
                    form16 = "log(totCover) ~ log(lag.cover)+TmeanSpr1+TmeanSpr2",
                    form17 = "log(totCover) ~ log(lag.cover)+ppt1*TmeanSpr1",
                    form18 = "log(totCover) ~ log(lag.cover)+ppt2*TmeanSpr2",
                    form19 = "log(totCover) ~ log(lag.cover)+ppt1+TmeanSpr1",
                    form20 = "log(totCover) ~ log(lag.cover)+ppt2+TmeanSpr2",
                    form21 = "log(totCover) ~ log(lag.cover)+ppt1",
                    form22 = "log(totCover) ~ log(lag.cover)+TmeanSpr1",
                    form23 = "log(totCover) ~ log(lag.cover)+ppt2",
                    form24 = "log(totCover) ~ log(lag.cover)+TmeanSpr2")


####
#### PRELIMINARIES --------------------------------
####

#load libraries
library(INLA)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../weather/Climate.csv")

#rename columns as necessary
# colnames(allD)[1] <- "year"

####
#### LOOP SPECIES ---------------------------------
####
dics <- matrix(ncol=length(sppList), nrow=length(fixed.forms))
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover  variable
  tmp=sppD[,c("quad","year","totCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  
  #Loop through models
  #set up unchanging random effects
  random <- "+f(quad, model='iid', prior='normal',param=c(0,0.001))+
    f(year, model='iid', prior='normal',param=c(0,0.001))"
  
  for(j in 1:length(fixed.forms)){
    fixed <- fixed.forms[[j]]
    formula1 <- as.formula(paste(fixed, random)) 
    out1 <- inla(formula1, data=growD,
                 family=c("gaussian"), verbose=TRUE,
                 control.compute=list(dic=T,mlik=T),
                 control.inla = list(h = 1e-10))
    dics[j,spp] <- out1$dic$dic
  }#end model loop
  
}#end species loop
  
#format dics and export
colnames(dics) <- sppList
dics <- as.data.frame(dics)
dics$modelNum <- c(1:length(fixed.forms))

outfile <- "./growthDIC.csv"
write.csv(dics, outfile)

 