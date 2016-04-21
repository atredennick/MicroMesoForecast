##  R Script to Process and Transform Raw Data into Dataframes Used
##  for Statistical Modeling and Model Simulations.
##
##  Author: Andrew Tredennick
##  Last update: 4-21-2016

# Clear the workspace
rm(list=ls())


####
####  LIBRARIES
####
library("plyr")
library("reshape2")



####
####  PRELIMINARIES
####
species <- c("BOGR", "HECO", "PASM", "POSE")
path2data <- "./speciesData/"
outpath <- "../processed_data/"
surv_outfile <- "surv_with_weather.RDS"
grow_outfile <- "grow_with_weather.RDS"
rec_outfile <- "rec_with_weather.RDS"
cover_outfile <- "cover_with_weather.RDS"
weather <- read.csv("../weather/Climate.csv")



####
####  SOURCE FILES
####
source("dataExtractQuad.R") # Adds in zeros for all measured quads; rms bad bogr quadyrs
source("format_quad_data_fxn.R") # Function for formatting quadrat data



####
####  QUADRAT PROPORTIONAL COVER DATA
####
quad_dat <- read.csv(paste0(path2data,"quadAllCover.csv"))
out_quad_dat <- format_quad_data(quad_dat, weather)
saveRDS(out_quad_dat, paste0(outpath,cover_outfile))



####
####  SURVIVAL DATA
####









