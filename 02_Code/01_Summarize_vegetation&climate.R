################################################################################
# Purpose: Summarize data from 51 shrub-grass interaction sites in sagebrush 
#          habitat across the western U.S.
#
#
# Rachel R. Renne
# March 26, 2025
# Updated: July 10, 2025
################################################################################

# Load libraries
library(dplyr)
library(daymetr)

# Set up directories
datadir <- '01_Data/Data for submission'
codedir <- '02_Code'
figdir <- '03_Figures/Figures for submission'

################################################################################
# Step 1: Get climate data for Sala et al. 1997 model

# Read in plots
plots <- read.csv(file.path(datadir, "plots.csv"))

# Remove climate variables if present to make data acquisition process clear
plots <- plots[,c(1:21)]

# Make day/month list
month <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),
           rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))

# Make dataframes:
ppt = data.frame(plot = plots$plot, Jan = NA, Feb = NA, Mar = NA, Apr = NA, 
                 May = NA, Jun = NA,Jul = NA, Aug = NA, Sep = NA, Oct = NA, 
                 Nov = NA, Dec = NA)
tmax = ppt
tmin = ppt

# Run through points and get data from daymet
for (i in 1:nrow(plots)){
  dmt <- download_daymet(site = "Daymet", lat = plots$latitude[i], 
                         lon = plots$longitude[i],start = 1991,end = 2020, 
                         internal = TRUE,silent = TRUE) 
  # Add month column
  dmt$data$month <- month 
  # Create data frame
  ppt[i,2:13] = tapply(dmt$data$prcp..mm.day., dmt$data$month, 
                       FUN = function(x){sum(x)/30})
  tmax[i,2:13] = tapply(dmt$data$tmax..deg.c., dmt$data$month, mean)
  tmin[i,2:13] = tapply(dmt$data$tmin..deg.c., dmt$data$month, mean)
}

# Put map and mat onto coords
plots$map <- apply(ppt[2:13], 1, sum)
plots$mat <- apply(cbind(tmin[,2:13],tmax[,2:13]), 1, mean)

# Check tmax and tmin are set up the same and calculate tavg
table(tmax$plot == tmin$plot)

# Set up tavg dataframe & calculate tavg from tmax and tmin
tavg <- tmax[0,]
for (i in 1:51){
  tavg[i,] <- rep(NA, 13)
  tavg$plot[i] <- tmax$plot[i]
  tavg[i,2:13] <- apply(rbind(tmax[i,2:13], tmin[i,2:13]), 2, mean)
}

# Change names to add to plot df
names(ppt)[2:13] <- paste0(names(ppt)[2:13],"_ppt")
names(tavg)[2:13] <- paste0(names(tavg)[2:13],"_tavg")

# Ensure plot order is the same
table(plots$plot == tavg$plot)
table(plots$plot == ppt$plot)

# Add onto plots df
plots1 <- merge(plots, ppt, by="plot")
plots <- merge(plots1, tavg, by="plot")

# Save updated plots DF
write.csv(plots, file.path(datadir,"plots.csv"), row.names = FALSE)

################################################################################
# Step 2: Summarize vegetation, soil water, and climate variables

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))

# Summarize climate
summary(plots$map)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   153.9   287.0   364.2   367.4   421.6   689.7

summary(plots$mat)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.917   5.863   7.464   7.340   9.273  10.670 

# Read in lpi data
lpi <- read.csv(file.path(datadir, "lpi.csv"))

# Create the ftype data frame
codes <- read.csv(file.path(datadir,"codes_withpathway.csv"))
codes1 <- codes[,c("code","GrowthForm")]
names(codes1) <- c("species","ftype")
ftype = codes1

# Calculate cover of each species
cover <- lpicover_ftypecorrected(lpi, ftype = ftype, n_layers = 5,
                  soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST","R"))

# Make dataframe wider
cover1 <- tidyr::pivot_wider(cover, id_cols = plot, names_from = species, 
                             values_from = cover, values_fill = 0)

# Look at shrub vs pgrass cover
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$PG ~ cover1$SH, pch = 16, col = rgb(0.5,0.5,0.5,0.7), cex = 2,
     xlab = "Shrub cover", ylab = "Perennial grass cover",
     xlim = c(0,0.72), ylim = c(0,0.72))

# Calculate shrub/grass dominance index:
cover1$dominance <- (cover1$SH - cover1$PG)/(cover1$SH + cover1$PG)

# Look at distribution of pft dominance
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
hist(cover1$dominance, xlim = c(-1,1), xlab = "Grass to Shrub dominance", 
     breaks = 10, main = "Functional type dominance")
abline(h = 0)
box()
abline(v = mean(cover1$dominance), lty = 2, lwd = 2, col = "red")

# Summarize cover
summary(cover1$SH)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09333 0.19667 0.25333 0.27451 0.33000 0.71333 

summary(cover1$PG)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0600  0.1267  0.1889  0.3067  0.5467 

# Summarize water balance
summary(plots$avg_trans_0_200cm)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 39.08  100.25  133.92  145.75  186.19  292.80 

summary(plots$prop_shallow)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3254  0.4355  0.4789  0.4847  0.5309  0.6909 

################################################################################
# Step 3: Summarize dung

# Read in dung
dung <- read.csv(file.path(datadir, "dung.csv"))

# Summarize by species per plot
dung1 <- dung %>% group_by(plot, animal) %>% summarise(total = sum(count))

# Calculate dung frequency per square meter
dung1$freq <- dung1$total/75

# Pivot
dung2 <- tidyr::pivot_wider(dung1, id_cols = plot, names_from = animal, 
                            values_from = freq, values_fill = 0)

# Summarize
summary(dung2$cow)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.01333 0.03946 0.06667 0.24000 

summary(dung2$horse)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.00000 0.00000 0.00000 0.01088 0.00000 0.14667 
