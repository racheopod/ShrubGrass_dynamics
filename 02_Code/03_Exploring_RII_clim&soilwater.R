################################################################################
# Purpose: Evaluate shrub-grass interactions across plots
#
#
# Rachel R. Renne
# April 2, 2025
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
# Step 1: Get cover variables

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))
source(file.path(codedir,"plantcommunity.R"))

# Read in plots file
plots <- read.csv(file.path(datadir, "plots.csv"))

# Read in lpi data
lpi <- read.csv(file.path(datadir, "lpi.csv"))

# Create the ftype data frame
codes <- read.csv(file.path(datadir,"codes_withpathway.csv"))
# Create dataframe convertain pathway codes to growthforms
pathdf <- data.frame(Code = unique(codes$Code),
                     GrowthForm = c("PG","PF","AF","SSH","SH","PG","AG","PG",
                                    "PG","SH","AG","NO GAPS","N","SUC","PG","PG",
                                    "PF","TR","AG","R","R","R","L","R","L","LC","M","PG",
                                    "S","R","VL","L"))
pathdf1 <- pathdf[!is.na(pathdf$GrowthForm),]

# Merge pathdf1 onto codes
codes1a <- merge(codes, pathdf1, by = "Code")
codes1 <- codes1a[,c("code","GrowthForm")]
names(codes1) <- c("species","ftype")
# Remove "S" (soil), "R" (rock), "N" (none), and "L" (litter)  rows
codes2 <- codes1[!(codes1$ftype %in% c("S","R","N","L")),]
ftype = codes2

# Calculate cover of each species
cover <- lpicover_ftypecorrected(lpi, ftype, n_layers = 5,
                                 soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"))

# Make dataframe wider
cover1 <- tidyr::pivot_wider(cover, id_cols = plot, names_from = species, 
                             values_from = cover, values_fill = 0)

# Look at some relationships
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$PG ~ cover1$SH, pch = 16, col = rgb(0.5,0.5,0.5,0.7), cex = 2)

# Calculate shrub/grass dominance index:
cover1$dominance <- (cover1$SH - cover1$PG)/(cover1$SH + cover1$PG)


################################################################################
# Step 2: Get RIIs for PG, AG, PF, AF

# Read in community data
pcomm <- read.csv(file.path(datadir, "community.csv"))

# Calculate both microsite and canopy condition values
ms_and_ca_calcs <- pcomm_microsite(pcomm, ftype, can_edge = FALSE)

# Separate out each one
microsites <- ms_and_ca_calcs[[1]]
canopy_condition <- ms_and_ca_calcs[[2]]

# Remove microsite records where there are 0 in canopy and microsite interspaces
#microsites1 <- microsites[microsites$canopy_n_ind > 0 & microsites$interspace_n_ind > 0, ]
# Or maybe don't remove these
microsites1 <- microsites

# Calculate some RIIs (cover, number)
microsites1$pcomm_cover_rii <- calculate_rii(canopy=microsites1$canopy_avg_cover, 
                                             interspace=microsites1$interspace_avg_cover)

microsites1$pcomm_number_rii <- calculate_rii(canopy=microsites1$canopy_avg_number, 
                                              interspace=microsites1$interspace_avg_number)

microsites1$pcomm_size_rii <- calculate_rii(canopy=microsites1$canopy_size,
                                            interspace=microsites1$interspace_size)

microsites1$pcomm_height_rii <- calculate_rii(canopy=microsites1$canopy_height,
                                              interspace=microsites1$interspace_height)
################################################################################
# Step 3: Get soil water variables

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

################################################################################
# Step 4: Calculate interactions

# Pull out individual ftypes for comparison
pgrass <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype == "PG",]
# Combine with soil water results
pgrass <- merge(pgrass, results, by = "plot", all = TRUE)

agrass <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype == "AG",]
# Combine with soil water results
agrass <- merge(agrass, results, by = "plot", all = TRUE)

pforb <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype == "PF",]
# Combine with soil water results
pforb <- merge(pforb, results, by = "plot", all = TRUE)

aforb <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype == "AF",]
# Combine with soil water results
aforb <- merge(aforb, results, by = "plot", all = TRUE)

aherb <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype %in% c("AG","AF"),]
pherb <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype %in% c("PG","PF"),]
allherb <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype %in% c("AG","AF","PG","PF"),]




# Summarize over plots for each
combineFT <- function(df){
  df <- df %>% 
    dplyr::group_by(plot) %>%
    dplyr::summarize(canopy_avg_number = sum(canopy_avg_number, na.rm = T),
                     interspace_avg_number = sum(interspace_avg_number, na.rm = T),
                     canopy_avg_cover = sum(canopy_avg_cover, na.rm = T),      
                     interspace_avg_cover = sum(interspace_avg_cover, na.rm = T),  
                     canopy_size = sum(canopy_size, na.rm = T),          
                     interspace_size = sum(interspace_size, na.rm = T),       
                     canopy_height = sum(canopy_height, na.rm = T),        
                     interspace_height = sum(interspace_height, na.rm = T))
  # Calculate some RIIs
  df$pcomm_cover_rii <- calculate_rii(canopy=df$canopy_avg_cover, 
                                      interspace=df$interspace_avg_cover)
  
  df$pcomm_number_rii <- calculate_rii(canopy=df$canopy_avg_number, 
                                       interspace=df$interspace_avg_number)
  
  df$pcomm_size_rii <- calculate_rii(canopy=df$canopy_size,
                                     interspace=df$interspace_size)
  
  df$pcomm_height_rii <- calculate_rii(canopy=df$canopy_height,
                                       interspace=df$interspace_height)
  return(df)
} 

# Combine fts
aherb <- combineFT(aherb)
# Combine with soil water results
aherb <- merge(aherb, results, by = "plot", all = TRUE)

pherb <- combineFT(pherb)
# Combine with soil water results
pherb <- merge(pherb, results, by = "plot", all = TRUE)

allherb <- combineFT(allherb)
# Combine with soil water results
allherb <- merge(allherb, results, by = "plot", all = TRUE)

################################################################################
# Step 4: Check for dung effects

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

# Add in plots with no dung
pt1 <- dung2[1:2,]
pt1$plot <- c("RkyMtn4-999", "GB5-005")
pt1[,2:8] <- 0
dung2a <- rbind(dung2, pt1)

# Merge with pgrass to check
dung3 <- merge(dung2a, pgrass, by = "plot", all = TRUE)

# Visualize 
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(dung3$pcomm_cover_rii~dung3$prop_shallow, pch = 16, 
     col = rgb(0.5,0.5,0.5,0.5),
     cex = 20*(dung3$cow+0.05),
     xlim = c(0.3,0.7), ylim = c(-1,1))

cor.test(dung3$pcomm_cover_rii, dung3$cow)
# r = -0.30, p = 0.034
cor.test(dung3$pcomm_cover_rii, dung3$horse)
# r = 0.101, p = 0.491
cor.test(dung3$pcomm_size_rii, dung3$cow)
# r = -0.257, p = 0.088
cor.test(dung3$pcomm_size_rii, dung3$horse)
# r = 0.186, p = 0.222
cor.test(dung3$pcomm_height_rii, dung3$cow)
# r = -0.165, p = 0.278
cor.test(dung3$pcomm_height_rii, dung3$horse)
# r = 0.149, p = 0.328
cor.test(dung3$pcomm_number_rii, dung3$cow)
# r = -0.227, p = 0.117
cor.test(dung3$pcomm_number_rii, dung3$horse)
# r = -0.045, p = 0.759

# Look at combined dung (none are significant)
cor.test(dung3$pcomm_cover_rii, apply(cbind(dung3$cow, dung3$horse), 1, sum))
cor.test(dung3$pcomm_size_rii, apply(cbind(dung3$cow, dung3$horse), 1, sum))
cor.test(dung3$pcomm_height_rii, apply(cbind(dung3$cow, dung3$horse), 1, sum))

# Use stepwise regression with prop_shallow & totaltransp & dung to predict pgrass interaction
mnull <- lm(dung3$pcomm_cover_rii ~ 1)
mfull <- lm(dung3$pcomm_cover_rii ~ dung3$prop_shallow*dung3$avg_trans_0_200cm*dung3$cow*dung3$horse)
msw_cow<- step(mnull, scope=list(upper=mfull),direction = "both",
               trace = 1)
summary(msw_cow)
par(mfrow = c(2,2), mgp = c(1,0.1,0), tcl = 0.1, mar = c(2,2,2,1))
plot(msw_cow)
# R2 = 0.4184, p = 1.833e-05


################################################################################
# Step 5: Create figures relating climate and soilwater to RIIs

# Merge pgrass with plots to get MAP and MAT
pgrass <- merge(pgrass, plots, by = "plot", all = T)
par(mfrow = c(1,1))

# A) MAP

# COVER
png(file.path(figdir, "Pgrass_RII_cover_vs_map.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_cover_rii, pgrass$map)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_cover_rii~pgrass$map, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual precipitation (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_cover_rii~pgrass$map, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# SIZE
png(file.path(figdir, "Pgrass_RII_size_vs_map.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_size_rii, pgrass$map)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_size_rii~pgrass$map, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual precipitation (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_size_rii~pgrass$map, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Number
png(file.path(figdir, "Pgrass_RII_number_vs_map.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_number_rii, pgrass$map)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_number_rii~pgrass$map, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual precipitation (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_number_rii~pgrass$map, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Height
png(file.path(figdir, "Pgrass_RII_height_vs_map.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_height_rii, pgrass$map)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_height_rii~pgrass$map, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual precipitation (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_height_rii~pgrass$map, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()



# B)  Total transpiration

# COVER
png(file.path(figdir, "Pgrass_RII_cover_vs_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_cover_rii, pgrass$avg_trans_0_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_cover_rii~pgrass$avg_trans_0_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_cover_rii~pgrass$avg_trans_0_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# SIZE
png(file.path(figdir, "Pgrass_RII_size_vs_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_size_rii, pgrass$avg_trans_0_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_size_rii~pgrass$avg_trans_0_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_size_rii~pgrass$avg_trans_0_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Number
png(file.path(figdir, "Pgrass_RII_number_vs_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_number_rii, pgrass$avg_trans_0_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_number_rii~pgrass$avg_trans_0_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_number_rii~pgrass$avg_trans_0_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Height
png(file.path(figdir, "Pgrass_RII_height_vs_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_height_rii, pgrass$avg_trans_0_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_height_rii~pgrass$avg_trans_0_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Mean annual transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_height_rii~pgrass$avg_trans_0_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()


# C)  Shallow proportion of transpiration

# COVER
png(file.path(figdir, "Pgrass_RII_cover_vs_prop_shallow.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_cover_rii, pgrass$prop_shallow)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_cover_rii~pgrass$prop_shallow, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "", cex = 1.4, ylim = c(-1,1))
mtext("Proportion of total transpiration\nfrom shallow soils (0-30 cm)",1,line = 2)
points(pgrass$pcomm_cover_rii~pgrass$prop_shallow, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p < 0.001"),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# SIZE
png(file.path(figdir, "Pgrass_RII_size_vs_prop_shallow.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_size_rii, pgrass$prop_shallow)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_size_rii~pgrass$prop_shallow, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "", cex = 1.4, ylim = c(-1,1))
mtext("Proportion of total transpiration\nfrom shallow soils (0-30 cm)",1,line = 2)
points(pgrass$pcomm_size_rii~pgrass$prop_shallow, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Number
png(file.path(figdir, "Pgrass_RII_number_vs_prop_shallow.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_number_rii, pgrass$prop_shallow)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_number_rii~pgrass$prop_shallow, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_number_rii~pgrass$prop_shallow, lwd = 1.7, col = 'black', cex = 1.4)
mtext("Proportion of total transpiration\nfrom shallow soils (0-30 cm)",1,line = 2)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Height
png(file.path(figdir, "Pgrass_RII_height_vs_prop_shallow.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_height_rii, pgrass$prop_shallow)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_height_rii~pgrass$prop_shallow, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "", cex = 1.4, ylim = c(-1,1))
mtext("Proportion of total transpiration\nfrom shallow soils (0-30 cm)",1,line = 2)
points(pgrass$pcomm_height_rii~pgrass$prop_shallow, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()




# D) Total shallow soil transpiration

# COVER
png(file.path(figdir, "Pgrass_RII_cover_vs_shallowtranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_cover_rii, pgrass$avg_trans_0_30cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_cover_rii~pgrass$avg_trans_0_30cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average shallow transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_cover_rii~pgrass$avg_trans_0_30cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# SIZE
png(file.path(figdir, "Pgrass_RII_size_vs_shallowtranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_size_rii, pgrass$avg_trans_0_30cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_size_rii~pgrass$avg_trans_0_30cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average shallow transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_size_rii~pgrass$avg_trans_0_30cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Number
png(file.path(figdir, "Pgrass_RII_number_vs_shallowtranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_number_rii, pgrass$avg_trans_0_30cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_number_rii~pgrass$avg_trans_0_30cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average shallow transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_number_rii~pgrass$avg_trans_0_30cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Height
png(file.path(figdir, "Pgrass_RII_height_vs_shallowtranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_height_rii, pgrass$avg_trans_0_30cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_height_rii~pgrass$avg_trans_0_30cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average shallow transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_height_rii~pgrass$avg_trans_0_30cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# E) Total deep soil transpiration

# COVER
png(file.path(figdir, "Pgrass_RII_cover_vs_deeptranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_cover_rii, pgrass$avg_trans_30_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_cover_rii~pgrass$avg_trans_30_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average deep transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_cover_rii~pgrass$avg_trans_30_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# SIZE
png(file.path(figdir, "Pgrass_RII_size_vs_deeptranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_size_rii, pgrass$avg_trans_30_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_size_rii~pgrass$avg_trans_30_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average deep transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_size_rii~pgrass$avg_trans_30_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Number
png(file.path(figdir, "Pgrass_RII_number_vs_deeptranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_number_rii, pgrass$avg_trans_30_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_number_rii~pgrass$avg_trans_30_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average deep transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_number_rii~pgrass$avg_trans_30_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# Height
png(file.path(figdir, "Pgrass_RII_height_vs_deeptranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(pgrass$pcomm_height_rii, pgrass$avg_trans_30_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(pgrass$pcomm_height_rii~pgrass$avg_trans_30_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "RII: Competitive (-) to facilitative (+)", 
     xlab = "Average deep transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(pgrass$pcomm_height_rii~pgrass$avg_trans_30_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

cor.test(pgrass$pcomm_cover_rii, pgrass$mat)
# r = 0.31, p = 0.030

################################################################################
# Step 5b: Look at correlation of RIIs for other plant types with each variable

checkriis <- function(x, y, ftype){
  print(ftype)
  for (ii in 1:ncol(x)){
    print(names(x)[ii])
    for (i in 1:ncol(y)){
      print(paste0("          ",names(y)[i]))
      mcor <- cor.test(x[,ii], y[,i])
      print(paste0("              r = ",round(mcor$estimate,2),", p = ", round(mcor$p.value, 5)))
    }
  }
}

# Pgrass
checkriis(pgrass[,17:20], pgrass[,c(21:24,37,38)], ftype = "Perennial grass")
 #Pgrass cover only
checkriis(pgrass[,c(17,17)], pgrass[,c(21:24,37,38)], ftype = "Perennial grass")


#  Pforb
checkriis(pforb[,17:20], pforb[,21:24], ftype = "Perennial forb")


#  Aforb
checkriis(aforb[,17:20], aforb[,21:24], ftype = "Annual forb")


# Agrass
checkriis(agrass[,17:20], agrass[,21:24], ftype = "Annual grass")


#  Aherb
checkriis(aherb[,10:13], aherb[,14:17], ftype = "Annual herbaceous")

plot(aherb$pcomm_height_rii~aherb$prop_shallow, 
     xlim = c(0.3,0.7), ylim = c(-1,1))
abline(h = 0, lty = 2)

# Compare height interaction pgrass and aherb
plot(aherb$pcomm_height_rii~pgrass$pcomm_height_rii,
     xlim = c(-1,1), ylim = c(-1,1))
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)


#  Pherb
checkriis(pherb[,10:13], pherb[,14:17], ftype = "Perennial herbaceous")


#  All herb
checkriis(allherb[,10:13], allherb[,14:17], ftype = "All herbaceous")


################################################################################
# Side note: Check on "forb" data for possible shrub-grass interactions
forbs <- read.csv(file.path(datadir, "forbs.csv"))

# DryGB2-019: C = 26, I = 4, RII_number = 0.75 (prop_shallow = 0.47)
# DryGB3-999: C = 11, I = 7, RII_number = 0.22 (prop_shallow = 0.43)
# GB5-005: C = 15, I = 4, RII_number = 0.58 (prop_shallow = 0.32)

# CP3-998
# GB3-014
# GB4-020

# Pretty much as expected, this is promising
# Need to systematically rate these images, enter data, and complete these calculations

################################################################################
# Step 6: Final figure showing RII_cover ~ amount + depth distribution of soil water

plot(pgrass$pcomm_cover_rii~pgrass$pcomm_number_rii, 
     xlim = c(-1,1), ylim = c(-1,1))
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(0,1,lty = 3)

plot(pgrass$pcomm_cover_rii~pgrass$pcomm_size_rii, 
     xlim = c(-1,1), ylim = c(-1,1))
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(0,1,lty = 3)

plot(pgrass$pcomm_cover_rii~pgrass$pcomm_height_rii, 
     xlim = c(-1,1), ylim = c(-1,1))


# Create binary variable for wet vs dry sites
mean(pgrass$avg_trans_0_200cm)
pgrass$dry <- as.numeric(pgrass$avg_trans_0_200cm < mean(pgrass$avg_trans_0_200cm))

# Create a second binary variable for wet vs. dry sites based on median
pgrass$dry_median <- as.numeric(pgrass$avg_trans_0_200cm < median(pgrass$avg_trans_0_200cm))

# Now create variable with several levels (1 = very dry, 2 = moderate, 3 = very wet)
pgrass$moisture <- NA
pgrass[pgrass$avg_trans_0_200cm < quantile(pgrass$avg_trans_0_200cm, 0.33),]$moisture <- 1
pgrass[pgrass$avg_trans_0_200cm >= quantile(pgrass$avg_trans_0_200cm, 0.66),]$moisture <- 2
pgrass[pgrass$avg_trans_0_200cm >= quantile(pgrass$avg_trans_0_200cm, 0.33) &
          pgrass$avg_trans_0_200cm < quantile(pgrass$avg_trans_0_200cm, 0.66)  ,]$moisture <- 3
# Set as a factor
pgrass$moisture <- as.factor(pgrass$moisture)

summary.aov(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$moisture))
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$moisture))

summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$dry))

summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$dry_median))

summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$avg_trans_0_200cm))

summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+as.factor(pgrass$dry)))

# Look at 3 moisure levels
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(pgrass$pcomm_cover_rii~pgrass$prop_shallow, pch = 16, 
     col = c(rgb(1,0.1,0,0.5), rgb(0.5,0.5,0.5,0.5), rgb(0,0.1,1,0.5))[pgrass$moisture],
     cex = 20*(dung2$cow+0.05),
     xlim = c(0.3,0.7), ylim = c(-1,1))
abline(h=0, lty = 2)
moist1 <- smooth.spline(x = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$moisture == 1, ]$prop_shallow, 
                        y = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$moisture == 1, ]$pcomm_cover_rii,
                        all.knots = FALSE, nknots = 10)
lines(moist1, col = "darkred", lwd=2)
moist2 <- smooth.spline(x = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$moisture == 2, ]$prop_shallow, 
                        y = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$moisture == 2, ]$pcomm_cover_rii,
                        all.knots = FALSE, nknots = 8)
lines(moist2, col = "darkgrey", lwd=2)
moist3 <- smooth.spline(x = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$moisture == 3, ]$prop_shallow, 
                        y = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$moisture == 3, ]$pcomm_cover_rii,
                        all.knots = FALSE, nknots = 10)
lines(moist3, col = "darkblue", lwd=2)


# Look at just wet and dry plots
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(pgrass$pcomm_cover_rii~pgrass$prop_shallow, pch = 16, 
     col = c(rgb(0,0.1,1,0.7), rgb(1,0.1,0,0.7))[pgrass$dry+1],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(-1,1))
abline(h=0, lty = 2)
dry1 <- smooth.spline(x = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$dry == 1, ]$prop_shallow, 
                      y = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$dry == 1, ]$pcomm_cover_rii,
                      all.knots = FALSE, nknots = 10, spar = 0.7)
lines(dry1, col = "darkred", lwd=2)
dry2 <- smooth.spline(x = pgrass[!is.na(pgrass$pcomm_cover_rii) & !is.na(pgrass$pcomm_cover_rii) & pgrass$dry == 0, ]$prop_shallow, 
                      y = pgrass[!is.na(pgrass$pcomm_cover_rii) & pgrass$dry == 0, ]$pcomm_cover_rii,
                      all.knots = FALSE, nknots = 10, spar = 0.7)
lines(dry2, col = "darkblue", lwd=2)

# Now with linear models only
m <- lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$dry)

png(file.path(figdir,"Chapter1_Fig5b.png"), width = 4, height = 4, units = "in", res = 300)
par(mar= c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(pgrass$pcomm_cover_rii~pgrass$prop_shallow, pch = 16, 
     col = c("#01665eb3", "#8c510ab3")[pgrass$dry+1],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(-1,1),
     ylab =  "RII: Competitive (-) to Facilitative (+)", 
     xlab = "")
#points(cover1$dominance~results$prop_shallow, cex = 2, lwd = 2)
mtext("Proportion of total transpiration from\nshallow soils (0-30 cm)", 1, line = 2)
abline(h=0, lty = 2)
curve(m$coefficients[1] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#01665e") # wet plots
curve(m$coefficients[1] + m$coefficients[3] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#8c510a") # dry plots
mtext("             Dry: Y = 1.91 - 3.6*X", side = 1, line = -2.2, adj = 0)
mtext("             Wet: Y = 1.67 - 3.6*X", side = 1, line = -1, adj = 0)
# Wet
points(x = 0.325, y = -1.009, cex = 2, col = "#01665eb3", pch = 16)
#points(x = 0.325, y = -1.009, cex = 2, lwd = 2)
lines(x = c(0.3,0.35), y = rep(-1.009,2), lwd = 3, col = "#01665e")
# Dry
points(x = 0.325, y = -0.86, cex = 2, col = "#8c510ab3", pch = 16)
#points(x = 0.325, y = -0.86, cex = 2, lwd = 2)
lines(x = c(0.3,0.35), y = rep(-0.86,2), lwd = 3, col = "#8c510a")
text(expression(R^2 ~ "= 0.32"), x = 0.67, y = 0.97)
text("p < 0.001", x = 0.67, y = 0.84)
dev.off()

png(file.path(figdir,"Chapter1_Fig5a.png"), width = 4, height = 4, units = "in", res = 300)
par(mar= c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(pgrass$pcomm_cover_rii~pgrass$map, pch = 16, 
     col = rgb(0.5,0.5,0.5,0.5),
     cex = 2,
     ylim = c(-1,1),
     ylab =  "RII: Competitive (-) to Facilitative (+)", 
     xlab = "")
points(pgrass$pcomm_cover_rii~pgrass$map, cex = 2, lwd = 2)
mtext("Mean annual precipitation (mm)", 1, line = 1.5)
abline(h=0, lty = 2)
mtext("  r = -0.2, p = 0.174", side = 1, line = -1.2, adj = 0)
dev.off()






#######################
# Try alternative interaction metric based on "forbs" dataset from 2022?