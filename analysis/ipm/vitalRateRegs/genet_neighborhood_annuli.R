##  This script uses RGEOS functions to calculate the area of
##    neighboring species within annuli of each focal genet in
##    the genet_xy input file. This replaces earlier versions of
##    species crowding to model species interactions. The program
##    loops through a focal species list to estimate crowding
##    and outputs a .csv file with the results. The results are then
##    used in vital rate regressions for survival and growth.

##  Author:       Peter Adler (updated by Andrew Tredennick)
##  Email:        atredenn@gmail.com
##  Last update:  4.27.2015 -- Updated Peter's version to work with Montana
##                             data.


## Clear the workspace
rm(list=ls(all=TRUE))

####
####  Load libraries ----------------------------------------------------------
####
library(rgeos)
library(maptools)
library(sp)


####
#### Define program functions -------------------------------------------------
####
makeDiscs <- function(x,y,annuli){
  centroid<-readWKT(paste("POINT(",x," ",y,")",sep=""))
  n <- length(annuli)
  output <- list(n)
  for(ii in 1:n){
    output[[ii]] <-gBuffer(centroid,width=annuli[ii],quadsegs=50)
  }
  return(output)
}

ringXquad <- function(discList){
  output <- numeric(length(discList))
  qpoly <- readWKT("POLYGON((0 0,100 0,100 100,0 100,0 0))")
  for(ii in 1:length(discList)){
    output[ii]<-gArea(gIntersection(qpoly,discList[[ii]]))/gArea(discList[[ii]])
  }
  return(output) 
}


####
####  Begin species looping and file set ups ----------------------------------
####
tmp_files <- list.files("../../speciesData")
tmp_files <- tmp_files[-grep("*.csv", tmp_files)]
species_list <- tmp_files[2:4]

shapefile_path <- "/Volumes/adlerlab/group/montanachart/lifetables/polys/Species/"
species_shapefiles <- c(paste(shapefile_path,"Bouteloua gracilis/", sep=""),
                        paste(shapefile_path,"Hesperostipa comata/", sep=""),
                        paste(shapefile_path,"Pascopyrum smithii/", sep=""),
                        paste(shapefile_path,"Poa secunda/", sep=""))
nbSpp <- species_list  # make sure these are listed in same order as sppShapes

outfile_path <- "../../speciesData/"

path_to_spatial_data <- "/Volumes/adlerlab/group/montanachart/ipm/speciesData/"

annuli <- c(seq(2,20,2),seq(25,50,5),seq(60,150,10))
Nrings <- length(annuli)

for(do_species in species_list){
  focalSppName <- do_species
  spatial_file <- paste("/",do_species, "_genet_xy.csv", sep="")
  focalSpp <- paste(path_to_spatial_data, do_species, spatial_file, sep="")
  outfile_name <- paste(do_species, "/", do_species, "_nbhood_rings.csv", sep="")
  outfile <- paste(outfile_path, outfile_name, sep="")
  
  # set up data frame for output
  template1 <- data.frame(species=character(0),quad=character(0),year=numeric(0),
                          genetID=numeric(0),stringsAsFactors=F)
  template2 <- as.data.frame(matrix(0,0,(length(nbSpp)+1)*Nrings))
  names(template2) <- paste(c(sort(rep(nbSpp,Nrings)),rep("propQ",Nrings)),c(0,annuli[1:(length(annuli)-1)]),annuli,sep=".")
  outD <-cbind(template1,template2)
  
  # import focal genets
  focalD <- read.csv(focalSpp)
  
  # list unique quad x year combinations
  quadyear <- unique(focalD[,c("quad","year")],MARGIN=2)
  
  # loop through quads and years
  for(i in 1:NROW(quadyear)){
    doQuad=quadyear$quad[i] ; doYear=quadyear$year[i]
    tmpD <- subset(focalD,  quad==doQuad & year==doYear)
    
    # import shapefiles
    polys <- list(length(nbSpp)) 
    for(j in 1:length(species_shapefiles)){
      tmpDir<-paste(species_shapefiles[j],doQuad,"/",sep="" )
      tmpfile<-paste(doQuad,"_",doYear,".shp",sep="" )
      if(sum(is.element(list.files(tmpDir),tmpfile))==1){
        polys[[j]]<-readShapePoly(paste(tmpDir,tmpfile,sep=""))
      }else{
        polys[[j]]<-0
      }
    } # next j
    
    
    # loop through focal genets
    for(k in 1:NROW(tmpD)){
      
      ID <- tmpD$trackID[k]; x <- tmpD$x[k] ; y <- tmpD$y[k] ; 
      
      # format output
      tout1 <- template1
      tout1[1,1:2] <- c(focalSppName,as.character(doQuad)); tout1[1,3:4] <- c(doYear,ID)
      tout2 <- template2; tout2[1,]<-0
      
      # make circles
      discList <- makeDiscs(x,y,annuli)
      
      # calculate proportion of each ring inside the quadrat
      propRing <- ringXquad(discList)
      
      # get areas
      nbarea <- numeric(Nrings)
      for(j in 1:length(species_shapefiles)){ 
        if(is.numeric(polys[[j]])==F){
          for(m in 1:Nrings){
            A <- gIntersection(discList[[m]],polys[[j]])
            if(!is.null(A)){
              nbarea[m] <- gArea(A) # calculate area only if A is not NULL
            }else{
              nbarea[m] <- 0
            }
          } # next m  
          if(nbSpp[j]==focalSppName){
            nbarea <- nbarea - tmpD$area[k] # subtract focal genet
            nbarea[nbarea<0] <- 0  # fix rounding error
          }  # end if 
          nbarea <- nbarea-c(0,nbarea[1:(Nrings-1)]) # calculate ring-specific areas  
          tout2[1,(1+(j-1)*Nrings):(Nrings+(j-1)*Nrings)] <- nbarea  
        } # end if polys null
      } # next j  
      
      tout2[1,(NCOL(tout2)-Nrings+1):NCOL(tout2)] <- propRing
      tout1 <- cbind(tout1,tout2)
      outD <- rbind(outD,tout1)
      
    } # next k  
    
    print(paste(i,"of",NROW(quadyear),"quad years complete",sep=" "))
    flush.console()
    
  } # next i
  
  write.csv(outD,outfile,row.names=F)
  
} # next species

