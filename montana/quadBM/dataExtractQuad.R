# Script to get data ready for quad-based model

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#


#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(reshape2)
library(plyr)

sppList <- list.files("../speciesData/")[1:4]

allRecs <- list()
for(i in 1:length(sppList)){
  doSpp <- sppList[i]
  tmpF <- paste("../speciesData/",doSpp,"/quadratCover.csv", sep="")
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

#write the file
outfile <- "../speciesData/quadAllCover2.csv"
write.csv(finalD, outfile)


