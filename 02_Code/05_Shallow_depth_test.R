################################
# Purpose: Shallow soil depth sensitivity analysis for RII ~ prop_shallow results
# 
# Rachel R. Renne
# February 19, 2025
##################################

# Load relevant libraries
library("rSFSW2")
library("rSOILWAT2")
library("tidyr")
library(viridis)

# Set up directories
datdir <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data"
figdir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/03_Figures/Final_figures'
codedir <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/02_Code"

################################################################################
# Step 1: Get SOILWAT output from interspace simulations

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
plotlist1 <- read.csv(file.path(datdir, "plots.csv"))
plotlist <- plotlist1$plot

# Make results df
results <- data.frame(plot = plotlist,avg_trans_0_10cm = NA,avg_trans_0_20cm = NA,
                      avg_trans_0_30cm = NA,avg_trans_0_40cm = NA,
                      avg_trans_0_50cm = NA,avg_trans_0_60cm = NA,
                      avg_trans_0_70cm = NA,avg_trans_0_80cm = NA,
                      avg_trans_0_90cm = NA,avg_trans_0_100cm = NA,
                      avg_trans_0_200cm = NA)

# Loop through and get relevant results
for (i in 1:nrow(results)){
  # Get for this plot
  thisplot <- trans[trans$Label == plotlist[i],52:417]
  
  # Get shallow, mid, deep
  results$avg_trans_0_200cm[i] <- sum(thisplot[1:15,])
  results$avg_trans_0_10cm[i] <- sum(thisplot[1:2,])
  results$avg_trans_0_20cm[i] <- sum(thisplot[1:3,])
  results$avg_trans_0_30cm[i] <- sum(thisplot[1:4,])
  results$avg_trans_0_40cm[i] <- sum(thisplot[1:5,])
  results$avg_trans_0_50cm[i] <- sum(thisplot[1:6,])
  results$avg_trans_0_60cm[i] <- sum(thisplot[1:7,])
  results$avg_trans_0_70cm[i] <- sum(thisplot[1:8,])
  results$avg_trans_0_80cm[i] <- sum(thisplot[1:9,])
  results$avg_trans_0_90cm[i] <- sum(thisplot[1:10,])
  results$avg_trans_0_100cm[i] <- sum(thisplot[1:11,])
}


# Calculate ratio of shallow to total transpiration for each definition of shallow
results$prop_shallow_10 <- results$avg_trans_0_10cm/results$avg_trans_0_200cm
results$prop_shallow_20 <- results$avg_trans_0_20cm/results$avg_trans_0_200cm
results$prop_shallow_30 <- results$avg_trans_0_30cm/results$avg_trans_0_200cm
results$prop_shallow_40 <- results$avg_trans_0_40cm/results$avg_trans_0_200cm
results$prop_shallow_50 <- results$avg_trans_0_50cm/results$avg_trans_0_200cm
results$prop_shallow_60 <- results$avg_trans_0_60cm/results$avg_trans_0_200cm
results$prop_shallow_70 <- results$avg_trans_0_70cm/results$avg_trans_0_200cm
results$prop_shallow_80 <- results$avg_trans_0_80cm/results$avg_trans_0_200cm
results$prop_shallow_90 <- results$avg_trans_0_90cm/results$avg_trans_0_200cm
results$prop_shallow_100 <- results$avg_trans_0_100cm/results$avg_trans_0_200cm

################################################################################
# Step 2: Calculate RIIs with plant community data 

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))
source(file.path(codedir, "plantcommunity.R"))

# Create the ftype data frame
codes <- read.csv(file.path(datdir,"codes_withpathway.csv"))
codes1 <- codes[,c("code","GrowthForm")]
names(codes1) <- c("species","ftype")
# Remove "S" (soil), "R" (rock), "N" (none), and "L" (litter)  rows
codes2 <- codes1[!(codes1$ftype %in% c("S","R","N")),]
ftype = codes2

# Read in community data
pcomm <- read.csv(file.path(datdir, "community.csv"))

# Calculate both microsite and canopy condition values
ms_and_ca_calcs <- pcomm_microsite(pcomm, ftype, can_edge = FALSE)

# Separate out each one (get just microsites)
microsites <- ms_and_ca_calcs[[1]]
# canopy_condition <- ms_and_ca_calcs[[2]]

# Remove microsite records where there are 0 in canopy and microsite interspaces
#microsites1 <- microsites[microsites$canopy_n_ind > 0 & microsites$interspace_n_ind > 0, ]
# Or maybe don't remove these
microsites1 <- microsites

# Calculate some RIIs
microsites1$pcomm_cover_rii <- calculate_rii(canopy=microsites1$canopy_avg_cover, 
                                             interspace=microsites1$interspace_avg_cover)

microsites1$pcomm_number_rii <- calculate_rii(canopy=microsites1$canopy_avg_number, 
                                              interspace=microsites1$interspace_avg_number)

microsites1$pcomm_size_rii <- calculate_rii(canopy=microsites1$canopy_size,
                                            interspace=microsites1$interspace_size)

microsites1$pcomm_height_rii <- calculate_rii(canopy=microsites1$canopy_height,
                                              interspace=microsites1$interspace_height)
names(microsites1)[2] <- "species"

# Separate by PFT
pg <- microsites1[microsites1$species == "PG" | is.na(microsites1$species),]
# pf <- microsites1[microsites1$species == "PF",]
# AF <- microsites1[microsites1$species == "AF",]
# AG <- microsites1[microsites1$species == "AG",]

# Merge with cover
results1 <- merge(results, pg[,c(1,17)], by = "plot", all.x = TRUE)

# Convert Nan's on Rii of cover to 0 (neutral--no pgrass anywhere)
results1[is.na(results1$pcomm_cover_rii),]$pcomm_cover_rii <- 0

# Make results figure
png(file.path(figdir,"Chapter1_Supp1.png"), width = 4, height = 4, units = "in", res = 300)
# Look at figure of RII ~ prop_shallow
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(results1$pcomm_cover_rii~results1$prop_shallow_10, 
     xlim = c(0, 1), ylim = c(-1,1),
     col = "white", pch = 16, cex = 2,
     xlab = "",
     ylab = "RII: Competition (-) to Facilitation (+)")
mtext("Proportion of total transpiration from\nshallow soils (0-30 cm)", 1, line = 2)
abline(h=0, lty = 2)

shalres <- data.frame(shallow = seq(10,100,by=10), coeff = NA, pvalue = NA)
for (i in 1:10){
  exp_results <- cor.test(results1$pcomm_cover_rii, results1[,12+i]) 
  shalres$coeff[i] <- exp_results$estimate
  shalres$pvalue[i] <- exp_results$p.value
  
  m <- lm(results1$pcomm_cover_rii~results1[,12+i])
  curve(m$coefficients[1]+x*m$coefficients[2], from = min(results1[,12+i]), to = max(results1[,12+i]), 
        col = viridis(10)[i], lwd = 2, add = TRUE)
}

legend("bottomleft", bty = "n", legend = seq(10,100,by=10), title = "Shallow depth (cm)",
       col = viridis(10), lwd = 2, cex = 0.7, y.intersp = 0.7)
legend("bottomright", bty = "n", legend = "p < 0.05 for all depths")
dev.off()

################################################################################
# Step 3: Calculate dominance and compare across different definitions of shallow 
#           soil water

# Read in lpi data
lpi <- read.csv(file.path(datdir, "lpi.csv"))


# Calculate cover of each species
cover <- lpicover_ftypecorrected(lpi, n_layers = 5,
                                 soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"))

# Make dataframe wider
cover1 <- tidyr::pivot_wider(cover, id_cols = plot, names_from = species, 
                             values_from = cover, values_fill = 0)

# Calculate shrub/grass dominance index:
cover1$dominance <- (cover1$SH - cover1$PG)/(cover1$SH + cover1$PG)

# Check for df alignment
table(cover1$plot == results$plot)

# Save figure of dominance with different definitions of shallow
png(file.path(figdir,"Chapter1_Supp2.png"), width = 4, height = 4, units = "in", res = 300)
# Look at figure of Dominance ~ prop_shallow
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~results1$prop_shallow_10, 
     xlim = c(0, 1), ylim = c(-1,1),
     col = "white", pch = 16, cex = 2,
     xlab = "",
     ylab = "Grass (-) to Shrub (+) dominance")
mtext("Proportion of total transpiration from\nshallow soils (0-30 cm)", 1, line = 2)
abline(h=0, lty = 2)

shalres <- data.frame(shallow = seq(10,100,by=10), coeff = NA, pvalue = NA)
for (i in 1:10){
  exp_results <- cor.test(cover1$dominance, results1[,12+i]) 
  shalres$coeff[i] <- exp_results$estimate
  shalres$pvalue[i] <- exp_results$p.value
  
  m <- lm(cover1$dominance~results1[,12+i])
  curve(m$coefficients[1]+x*m$coefficients[2], from = min(results1[,12+i]), to = max(results1[,12+i]), 
        col = viridis(10)[i], lwd = 2, add = TRUE, lty = c(2,1)[(shalres$pvalue[i] <= 0.05)+1])
}

legend("bottomleft", bty = "n", legend = seq(10,100,by=10), title = "Shallow depth (cm)",
       col = viridis(10), lwd = 2, cex = 0.7, y.intersp = 0.7)
legend("bottomright", bty = "n", legend = "p < 0.05 for all depths > 10 cm")
dev.off()
