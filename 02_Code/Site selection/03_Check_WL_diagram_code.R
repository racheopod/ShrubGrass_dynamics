#################################
#
# Purpose: Verify that each candidate sampling location has appropriate climate patterns
#
# Inputs: "/Renne2021_Points.csv"
#
# Outputs: Walter-Leith diagrams for all sample points
#
# Rachel R. Renne
# April 13, 2021
############################

library(daymetr)
library(climatol)

# Read in pts
pts <- read.csv("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Points/Renne2021_Points.csv", row.names = 1)

# Make day/month list
month <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))

# Run through points and make png of WL diagram:
for (i in 2:nrow(pts)){
dmt <- download_daymet(site = "Daymet", lat = pts$ycoord[i], lon = pts$xcoord[i],
                start = 1981,
                end = 2010,
                internal = TRUE,
                silent = FALSE) 

# Add month column
dmt$data$month <- month 

# Create dat
dat <- t(data.frame(ppt = tapply(dmt$data$prcp..mm.day., dmt$data$month, FUN = function(x){sum(x)/30}),
                    tmax = tapply(dmt$data$tmax..deg.c., dmt$data$month, mean),
                    tmax = tapply(dmt$data$tmin..deg.c., dmt$data$month, mean),
                    tabsmin = tapply(dmt$data$tmin..deg.c., dmt$data$month, min)))

# Set up file
png(paste0("C:/Users/Rachel/OneDrive - The University of Montana/Personal/WL Diagrams/", pts$PlotID[i],".png"), width = 6, height = 6, units = "in", res = 300)
# Draw WL
diagwl(dat, est = pts$PlotID[i], alt = dmt$altitude, per = "1981-2010", mlab = "en")
dev.off()
}       
