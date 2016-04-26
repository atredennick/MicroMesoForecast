##  R Script for Plotting Results from Climate Change Scenarios
##
##  Author: Andrew Tredennick
##  Last update: 4-26-2016

# Clear the workspace
rm(list=ls())



####
####  PRELIMINARIES
####
path2qbms <- "../analysis/quadBM/simulations/results/climchange_varyparams/"
path2ipms <- "../analysis/ipm/simulations/results/climchange_varyparams/"
path2figs <- "../manuscript/components/"



####
####  LOAD LIBRARIES
####
library("plyr") # data summarising
library("reshape2") # data wrangling
library("ggplot2") # plotting
library("ggthemes") # plotting themes