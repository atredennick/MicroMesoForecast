
##this script is to organize the climate data for Montana site, based on the raw data.

setwd("Y:\\montanachart\\records\\weather")
raw<-read.csv("daily_climate_data.csv")
##This data is from Lance, including 1931 records, a previous year of first mapping year (1932)

#####for precipitation
raw$precip2<-as.numeric(levels(raw$precip))[raw$precip]
raw$precip2[is.na(raw$precip2)]<-0.005 ##trace amounts of precipitation (Peter): 1077 records replaced
##remove original precipitation records
raw<-raw[,c("year","month","day","max","min","precip2")]
names(raw)<-c("year","month","day","max","min","precip")


######for temperature
##there are max and min temperature records here, so we firstly need to average them
raw$temp<-(raw$max+raw$min)/2 ##average
##remove max and min columns
raw<-raw[,c("year","month","day","temp","precip")]
         
    
#############################################
#############################################
###take separate treatments for temperature and precipitation
raw_temp<-raw[,c("year","month","day","temp")]
raw_precip<-raw[,c("year","month","day","precip")]


###aggregate temperature
month_temp<-aggregate(raw_temp[,c("year","month","temp")],by=list(raw_temp$month,raw_temp$year),mean)
month_temp<-month_temp[,c("year","month","temp")]

###aggregate precipitation
month_precip<-aggregate(raw_precip[,c("precip")],by=list(raw_precip$month,raw_precip$year),sum)
names(month_precip)<-c("month","year","precip")
month_precip<-month_precip[,c("year","month","precip")]

###then we merge temperature and precipitation together
month_climate<-merge(month_precip,month_temp,by=c("year","month"))
##order the data frame according to year and month
month_climate<-month_climate[order(month_climate$year,month_climate$month),]





###########################################
###########################################
###need to get the seanonal climate data: mainly on SPRING

###firstly work on individual metric, then merge them together with YEAR
##total annual precipitation (for pptLag)
precip_annual<-aggregate(month_climate$precip,by=list(month_climate$year),sum)
names(precip_annual)<-c("year","precip")

###spring temperature
spring<-c(4,5,6) #April, May and June
month_climate_spr<-subset(month_climate,month %in% spring) ##a subset data for SPRING
###summer temperature
summer<-c(6,7,8,9) #June, Junly, Aug. and Sept...
month_climate_sum<-subset(month_climate,month %in% summer) ##a subset data for SUMMER

###for spring temperature
temp_spr<-aggregate(month_climate_spr$temp,by=list(month_climate_spr$year),mean)
names(temp_spr)<-c("year","temp")
###for summer temperature
temp_sum<-aggregate(month_climate_sum$temp,by=list(month_climate_sum$year),mean)
names(temp_sum)<-c("year","temp")



##for spring precipitation
precip_spr<-aggregate(month_climate_spr$precip,by=list(month_climate_spr$year),sum)
names(precip_spr)<-c("year","precip")



###merge them together: one by one
seasonal_climate_spr<-merge(precip_spr,temp_spr,by=c("year"))
seasonal_climate_sum<-merge(precip_spr,temp_sum,by=c("year"))



names(seasonal_climate_spr)<-c("year","precip_spr","temp_spr")
names(seasonal_climate_sum)<-c("year","precip_spr","temp_sum")
seasonal_climate_spr_sum<-merge(seasonal_climate_spr,seasonal_climate_sum,by=c("year","precip_spr"))
names(seasonal_climate_spr_sum)<-c("year","precip","temp_spr","temp_sum")

seasonal_climate<-merge(precip_annual,seasonal_climate_spr_sum,by=c("year"))

names(seasonal_climate)<-c("year","precip_anuual","precip_spr","temp_spr","temp_sum")




###Then make the data frame of climatic variables of interest
##get pptLag
tmp_pptLag<-seasonal_climate[,c("year","precip_anuual")]
names(tmp_pptLag)<-c("year","pptLag")
tmp_pptLag$year<-tmp_pptLag$year+1

my_pptLag<-merge(tmp_pptLag,seasonal_climate)
my_pptLag<-my_pptLag[,c("year","pptLag")]
my_pptLag<-my_pptLag[-dim(my_pptLag)[1],]


##get ppt1 and ppt2
tmp_ppt12<-seasonal_climate[-1,c("year","precip_spr")] ##no need of 1931 record here
tmp_ppt1<-tmp_ppt12[1:(dim(tmp_ppt12)[1]-1),] ##precip of the first year of the transition
tmp_ppt2<-tmp_ppt12[2:dim(tmp_ppt12)[1],]     ##precip of the second year of the transition


##get TmeanSpr1 and TmeanSpr2--Spring
tmp_temp12<-seasonal_climate[-1,c("year","temp_spr")] ##no need of 1931 record here
tmp_temp1<-tmp_temp12[1:(dim(tmp_temp12)[1]-1),] ##temp of the first year of the transition
tmp_temp2<-tmp_temp12[2:dim(tmp_temp12)[1],]     ##temp of the second year of the transition

##get TmeanSum1 and TmeanSum2--Summer
tmp_temp34<-seasonal_climate[-1,c("year","temp_sum")] ##no need of 1931 record here
tmp_temp3<-tmp_temp34[1:(dim(tmp_temp34)[1]-1),] ##temp of the first year of the transition
tmp_temp4<-tmp_temp34[2:dim(tmp_temp34)[1],]     ##temp of the second year of the transition


###combine these variables together
##here the name 'Climate' follows Peter's
Climate<-my_pptLag
Climate$ppt1<-tmp_ppt1[,2]
Climate$TmeanSpr1<-tmp_temp1[,2]
Climate$ppt2<-tmp_ppt2[,2]
Climate$TmeanSpr2<-tmp_temp2[,2]
Climate$TmeanSum1<-tmp_temp3[,2]
Climate$TmeanSum2<-tmp_temp4[,2]

write.table(Climate,"Climate.csv",row.names=F,sep=",")



# ###for Revision JE: the data we used there is from 1915 to 1950 (here with the index of 1:36)
# png("Directional environment change.png",units="in",height=8,width=5,res=400)
# 
# par(mfrow=c(2,1))
# 
# plot(Climate$year,Climate$ppt1,xlab="Year",ylab="Mean annual precipitation")
# plot(Climate$year,Climate$TmeanSpr1,xlab="Year",ylab="Mean annual temperature")
# 
# dev.off()

