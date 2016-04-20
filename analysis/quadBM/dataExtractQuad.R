##  Sript to takes quadrat level data by species
##    and further cleans it for analysis by adding in zeros
##    for quads observed but not recorded, and taking out
##    "bad" BOGR quad-year

##  Loops through species-level data folders;
##    opens data and then adds in zeros for observed
##    quadrats that had no records due to zeros.

##  Date started: 1-5-2015
##  Date completed: 1-10-2015
##  Date updated: 3-20-2015 // New lines to take out bad BOGR quad-years

##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

####
####  Load libraries -------------------------------------
####
library(reshape2)
library(plyr)

####
#### Loop through species folders to get data -----------
####
#get species list
sppList <- list.files("../speciesData/")[1:4]

allRecs <- list()
for(i in 1:length(sppList)){
  doSpp <- sppList[i]
  if(doSpp=="BOGR"){ tmpF <- paste("../speciesData/",doSpp,"/edited/quadratCover.csv", sep="") }
  if(doSpp!="BOGR"){tmpF <- paste("../speciesData/",doSpp,"/quadratCover.csv", sep="") }
  tmpD <- read.csv(tmpF)
  tmpD$Species <- doSpp
  allRecs <- rbind(allRecs, tmpD)
}
allRecs$propCover <- allRecs$totCover/10000

#bring in quadrat inventory
quadInv <- read.csv("../speciesData/quad_inventory_2.csv")
quadInv$year <- quadInv$A1
quadM <- melt(quadInv, id.vars = "year")

#add in zeros for quads that were actually measured in a given year
# uQuads <- unique(allRecs$quad)
yrs <- unique(allRecs$year)
addLines <- data.frame(quad = NA, 
                       year = NA, 
                       totCover = NA, 
                       Species = NA, 
                       propCover = NA)

for (s in 1:length(sppList)){
  tmpDspp <- subset(allRecs, Species==sppList[s])
  doSpp <- sppList[s]
  for(y in 1:length(yrs)){
    tmpD <- subset(tmpDspp, year==yrs[y])
    tmpQ <- subset(quadM, year==yrs[y] & is.na(value)!=TRUE)
    uQuads <- tmpQ$variable
    for(q in 1:length(uQuads)){
      tmpQ <- uQuads[q]
      test <- which(tmpD$quad==as.character(tmpQ))
      ifelse(length(test)!=0,
             print(paste("quad ", tmpQ, " is good")),
             tmpLine <- data.frame(quad = uQuads[q], 
                                   year = yrs[y], 
                                   totCover = 0, 
                                   Species = doSpp, 
                                   propCover = 0))
      ifelse(exists("tmpLine")==TRUE,
             addLines <- rbind(addLines, tmpLine),
             print("nada"))
      rm(tmpLine)
    }
  }
}

#add in zeros
finalD <- rbind(allRecs, addLines)
finalD <- finalD[with(finalD, order(Species, year, quad)), ]
finalD <- finalD[1:(nrow(finalD)-1), ] #removes NA row

####
####  Remove bad BOGR quad-years -----------------------
####
# This reads from the Adler lab server to make sure we
#   get the same bad BOGR quad-years every time (didn't want to 
#   copy to this folder in case we update the original file).
#   These were taken out in the data prep stage, but added back in here
#   since we add in zeros based on the quadrat inventory.
bad_bogrs <- read.csv("/Volumes/adlerlab/group/montanachart/suspect_BOGR_quads.csv")
bad_bogrs$year <- bad_bogrs$year-1900 #to match finalD
bad_bogrs$bad_check <- "yes" #just a flaggin column
bad_bogr_data <- merge(finalD, bad_bogrs, by.x=c("quad", "year"), 
                       by.y = c("quadrat", "year"), all.x=TRUE)
bad_bogr_data <- bad_bogr_data[order(bad_bogr_data$Species, bad_bogr_data$year, bad_bogr_data$quad),]

#Make sure the ordering is right
sum(finalD$propCover-bad_bogr_data$propCover) #should return 0

#OK, now find rows to get rid of
ids_to_remove <- which(bad_bogr_data$bad_check == "yes" & bad_bogr_data$Species=="BOGR")
finalD <- finalD[-ids_to_remove,]

####
####  Write the file --------------------------------
####
outfile <- "../speciesData/quadAllCover.csv"
write.csv(finalD, outfile)

# plotdf <- ddply(finalD, .(quad,year,Species), summarise,
#                 meancov = mean(totCover))
# ggplot(plotdf, aes(year,meancov,color=Species))+
#   geom_line()+
#   facet_wrap("quad")
