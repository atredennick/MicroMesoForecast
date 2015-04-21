##  This script subsets out recruit area for quadrats that no longer
##    have recruit areas that can be fit because we removed the years
##    before that event. We took out quad-years for BOGR that looked
##    suspicious in terms of mapping. Now that we took those out, we
##    need to also update the recruitArea.csv file so that we do not
##    have recruit areas for what is now the first year of our observations.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4-21-2015

##  Note that this script requires a connection the Adler Lab USU server
##    for access to files that I don't want to make multiple copies of.

# Clear the workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(reshape2)
library(plyr)

####
####  Read in necessary data files ---------------------------------------------
####
bad_bogrs <- read.csv("/Volumes/adlerlab/group/montanachart/suspect_BOGR_quads.csv")
current_recruit_area <- read.csv("recArea_quads_not_removed.csv")
full_quad_inventory <- read.csv("/Volumes/adlerlab/group/montanachart/records/quad_inventory_2.csv")
suppressMessages(full_quad_long <- melt(full_quad_inventory))
colnames(full_quad_long) <- c("quad", "year")
rm(full_quad_inventory)


####
####  Find instances where initial year(s) were removed ------------------------
####
# Find minimum year for each quad
bad_bogrs_beginyear <- ddply(bad_bogrs, .(quadrat), summarise,
                             first_year = min(year))
inventory_beginyear <- ddply(full_quad_long, .(quad), summarise,
                             first_year = min(year, na.rm = TRUE))

