##  Sript to takes quadrat level data by species
##    and further cleans it for analysis by adding in zeros
##    for quads observed but not recorded

##  Loops through species-level data folders;
##    opens data and then adds in zeros for observed
##    quadrats that had no records due to zeros.

##  Date started: 1-5-2015
##  Date completed: 1-10-2015

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
sppList <- "PSSP"

allRecs <- list()
for(i in 1:length(sppList)){
  doSpp <- sppList[i]
  tmpF <- paste("./",doSpp,"/quadratCover.csv", sep="")
  tmpD <- read.csv(tmpF)
  tmpD$Species <- doSpp
  allRecs <- rbind(allRecs, tmpD)
}
allRecs$propCover <- allRecs$totCover/10000

#bring in quadrat inventory
quadInv <- read.csv("./quad_inventory.csv")
quadInv$year <- c(min(quadInv, na.rm = T):max(quadInv, na.rm = T))
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
      if(exists("tmpLine")==TRUE){
        addLines <- rbind(addLines, tmpLine)
        rm(tmpLine)
      } # end if
      
    }
  }
}

#add in zeros
finalD <- rbind(allRecs, addLines)
finalD <- finalD[with(finalD, order(Species, year, quad)), ]
finalD <- finalD[1:(nrow(finalD)-1), ] #removes NA row


####
####  Write the file --------------------------------
####
# outfile <- "./quadAllCover.csv"
# write.csv(finalD, outfile)


