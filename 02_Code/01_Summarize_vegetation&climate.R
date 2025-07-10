################################################################################
# Purpose: Summarize data from shrub-grass interaction plots
#
#
# Rachel R. Renne
# March 26, 2025
################################################################################

# Load libraries
library(dplyr)
library("rSFSW2")
library("rSOILWAT2")

# Set up directories
datadir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data'
codedir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/02_Code'
outdir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/01_Data'
figdir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/03_Figures/Final_figures'

################################################################################
# Step 1: Summarize vegetation and climate variables

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))

# Read in plots file
plots <- read.csv(file.path(datadir, "plots.csv"))

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
codes <- read.csv(file.path(datadir,"codes.csv"))
codes1 <- codes[,c("code","GrowthForm")]
names(codes1) <- c("species","ftype")
# Remove "S" (soil), "R" (rock), "N" (none), and "L" (litter)  rows
codes2 <- codes1[!(codes1$ftype %in% c("S","R","N","L")),]
ftype = codes2

# Calculate cover of each species
cover <- lpicover_ftypecorrected(lpi, n_layers = 5,
                  soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"))

# Make dataframe wider
cover1 <- tidyr::pivot_wider(cover, id_cols = plot, names_from = species, 
                             values_from = cover, values_fill = 0)

# Look at some relationships
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$PG ~ cover1$SH, pch = 16, col = rgb(0.5,0.5,0.5,0.7), cex = 2,
     xlab = "Shrub cover", ylab = "Perennial grass cover",
     xlim = c(0,0.72), ylim = c(0,0.72))

hist(cover1$PG, breaks = 15, main = "",
     xlab = "Perennial grass cover (%)",
     xlim = c(0,0.72))
abline(h=0)
box()

hist(cover1$SH, breaks = 15, main = "",
     xlab = "Shrub cover (%)",
     xlim = c(0,0.72))
abline(h=0)
box()


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

################################################################################
# Step 2: Summarize soil water variables

# Set up project directory
dir_prj <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/04_Simulation_Experiment/Interspace_simulations/Interspace_simulations_20250114"

# Database directory
dir_dat <- file.path(dir_prj, "4_Simulation")

# Database name
outDB_fname <- file.path(dir_dat, "dbOutput.sqlite3")

# Connect to output database
outDB_con <- RSQLite::dbConnect(RSQLite::SQLite(), outDB_fname)

# Average daily transpiration
trans <- get.Table_Scenario(outDB_fname,
                            responseName = "aggregation_doy_Transpiration", MeanOrSD = "Mean",
                            scenario = "Current", header = TRUE)

# Clean up
RSQLite::dbDisconnect(outDB_con)  

# Create a plotlist:
plotlist1 <- read.csv(file.path(dir_prj, "1_Input/SWRuns_InputMain_interspace_simulations_v12.csv"))
plotlist <- plotlist1$Label

# Make results df
results <- data.frame(plot = plotlist,avg_trans_0_30cm = NA,avg_trans_30_200cm = NA)

# Loop through and get relevant results
for (i in 1:nrow(results)){
  # Get for this plot
  thisplot <- trans[trans$Label == plotlist[i],52:417]
  
  # Get shallow, mid, deep
  results$avg_trans_0_30cm[i] <- sum(thisplot[1:4,])
  results$avg_trans_30_200cm[i] <- sum(thisplot[5:15,])
}

# Calculate total transpiration
results$avg_trans_0_200cm <- results$avg_trans_0_30cm+results$avg_trans_30_200cm

# Calculate ratio of shallow to total transpiration
results$prop_shallow <- results$avg_trans_0_30cm/results$avg_trans_0_200cm

# Summarize water balance
summary(results$avg_trans_0_200cm)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 39.08  100.25  133.92  145.75  186.19  292.80 

summary(results$prop_shallow)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3254  0.4355  0.4789  0.4847  0.5309  0.6909 

################################################################################
# Step 3: Summarize dung

# Read in dung
dung <- read.csv(file.path(datadir, "dung.csv"))

# Remove missing value
dung <- dung[!is.na(dung$count),]

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