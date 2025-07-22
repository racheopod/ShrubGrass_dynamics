################################
# Purpose: Shallow soil depth sensitivity analysis for RII ~ prop_shallow results
# 
# Rachel R. Renne
# February 19, 2025
# Updated: July 11, 2025
##################################

# Approach: We calculated the proportion of transpiration from shallow soils for
# shallow soils defined by depths 0 to: 10, 20, 30, 40, 50, 60, 70, 80, 90, and 100 cm.
# Here, we calculate the correlations between each set of prop_shallow results 
# and both RII and the dominance index.

# Load relevant libraries
library(tidyr)
library(viridis)

# Set up directories
datdir <- "01_Data/Data for submission"
codedir <- "02_Code"
figdir <- "03_Figures/Figures for submission"

################################################################################
# Step 1: Get shallow soil sensitivity analysis simulation results

results <- read.csv(file.path(datdir,"Shallow_soil_sensitivity_simulation_results.csv"))

################################################################################
# Step 2: Calculate RIIs with plant community data 

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))
source(file.path(codedir, "plantcommunity.R"))

# Create the ftype data frame
codes <- read.csv(file.path(datdir,"codes_withpathway.csv"))
codes1 <- codes[,c("code","GrowthForm")]
# Change GrowthForm for NO GAPS
codes1[is.na(codes1$GrowthForm),]$GrowthForm <- "NO GAPS"
names(codes1) <- c("species","ftype")
ftype = codes1

# Read in community data
pcomm <- read.csv(file.path(datdir, "community.csv"))

# Calculate both microsite and canopy condition values
ms_and_ca_calcs <- pcomm_microsite(pcomm, ftype, can_edge = FALSE)

# Pull out just microsites RII
# In this case, "canopy" microsite is canopy > 0, "interspace" is canopy = 0
microsites1 <- ms_and_ca_calcs[[1]]

# Calculate cover values for each microsite
microsites1$pcomm_cover_rii <- calculate_rii(canopy=microsites1$canopy_avg_cover, 
                                             interspace=microsites1$interspace_avg_cover)

# Separate out just perennial grasses
pg <- microsites1[microsites1$ftype == "PG" | is.na(microsites1$ftype),]

# Merge with cover
results1 <- merge(results, pg[,c(1,17)], by = "plot", all.x = TRUE)

# Make results figure
png(file.path(figdir,"Chapter1_Supp3b.png"), width = 4, height = 4, units = "in", res = 600)
# Look at figure of RII ~ prop_shallow
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(results1$pcomm_cover_rii~results1$prop_shallow_10, 
     xlim = c(0, 1), ylim = c(-1,1),
     col = "white", pch = 16, cex = 2,
     xlab = "",
     ylab = "RII: Competition (-) to Facilitation (+)")
mtext("Proportion of total transpiration from\nshallow soil layers (0-30 cm)", 1, line = 2)
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
cover <- lpicover_ftypecorrected(lpi, ftype, n_layers = 5,
                                 soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST","R"))

# Make dataframe wider
cover1 <- tidyr::pivot_wider(cover, id_cols = plot, names_from = species, 
                             values_from = cover, values_fill = 0)

# Calculate shrub/grass dominance index:
cover1$dominance <- (cover1$SH - cover1$PG)/(cover1$SH + cover1$PG)

# Check for df alignment
table(cover1$plot == results$plot)

# Save figure of dominance with different definitions of shallow
png(file.path(figdir,"Chapter1_Supp3a.png"), width = 4, height = 4, units = "in", res = 600)
# Look at figure of Dominance ~ prop_shallow
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~results1$prop_shallow_10, 
     xlim = c(0, 1), ylim = c(-1,1),
     col = "white", pch = 16, cex = 2,
     xlab = "",
     ylab = "Grass (-) to Shrub (+) dominance")
mtext("Proportion of total transpiration from\nshallow soil layers (0-30 cm)", 1, line = 2)
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
