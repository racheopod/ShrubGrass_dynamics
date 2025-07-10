################################################################################
# Purpose: Evaluate plant functional type dominance across plots
#
#
# Rachel R. Renne
# March 26, 2025
################################################################################

# Load libraries
library(dplyr)
library("rSFSW2")
library("rSOILWAT2")
library(terra)

# Set up directories
datadir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data'
codedir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/02_Code'
outdir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/01_Data'
figdir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/03_Figures/Final_figures'

################################################################################
# Step 1: Get cover variables

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))

# Read in plots file
plots <- read.csv(file.path(datadir, "plots.csv"))

# Read in lpi data
lpi <- read.csv(file.path(datadir, "lpi.csv"))

# Create the ftype data frame
codes <- read.csv(file.path(datadir,"codes_withpathway.csv"))
# Create dataframe convertain pathway codes to growthforms
pathdf <- data.frame(Code = unique(codes$Code),
                     GrowthForm = c("PG","PF","AF","SSH","SH","PG","AG","PG",
                                    "PG","SH","AG",NA,NA,"SUC","PG","PG",
                                    "PF","TR","AG",NA,NA,NA,NA,NA,NA,NA,NA,"PG",
                                    NA,NA,NA,NA))
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
# Step 2: Get climate data for Sala et al. 1997 model

# Read in precip and temp
ppt <- read.csv(file.path(datadir, "ppt.csv"))
tmax <- read.csv(file.path(datadir, "tmax.csv"))
tmin <- read.csv(file.path(datadir, "tmin.csv"))

# Calculate tavg
table(tmax$plot == tmin$plot)

# Set up tavg dataframe
tavg <- tmax[0,]
for (i in 1:51){
  tavg[i,] <- rep(NA, 13)
  tavg$plot[i] <- tmax$plot[i]
  tavg[i,2:13] <- apply(rbind(tmax[i,2:13], tmin[i,2:13]), 2, mean)
}

# Now calculate correlation
table(cover1$plot == tavg$plot & tavg$plot == ppt$plot)
cover1$CORR_PT <- NA
for (i in 1:51){
  m <- cor.test(t(ppt[i,2:13]), t(tavg[i,2:13]))
  cover1$CORR_PT[i] <- m$estimate
}

# Read in soil texture
st <- read.csv(file.path(datadir, "soiltexture.csv"))

# Calculate average texture of top 30 cm
cover1$sand <- NA
cover1$clay <- NA
for (i in 1:51){
  thissoil <- st[st$plot == cover1$plot[i] & st$depth_1 < 31 & !is.na(st$depth_1), ]
  cover1$sand[i] <- mean(thissoil$sand)
  cover1$clay[i] <- mean(thissoil$clay)
}

# Calculate water holding capacity (theta_S)
# theta_S = theta_33 + theta_S_33 - 0.097S + 0.043
# theta_33 = theta_33t + (1.283*(theta_33t^2) 0.374*theta_33t - 0.015)
# theta_33t = -0.251S + 0.195C + 0.011OM + 
#             0.006*(S*OM) - 0.027*(C*OM) +
#             0.452*(S*C) + 0.299
# theta_S_33 = theta_S_33t + (0.636*theta_S_33t - 0.107)
# theta_S_33t = 0.278*S + 0.034*C + 0.022*OM - 
#               0.018*(S*OM) - 0.027*(C*OM) -
#               0.584*(S*C) + 0.078


# Set OM to 0.0025
OM = 0.0025
theta_33t = -0.251*cover1$sand + 0.195*cover1$clay + 0.011*OM + 
             0.006*(cover1$sand*OM) - 0.027*(cover1$clay*OM) +
             0.452*(cover1$sand*cover1$clay) + 0.299

theta_33 = theta_33t + (1.283*(theta_33t^2) - 0.374*theta_33t - 0.015)

theta_S_33t = 0.278*cover1$sand + 0.034*cover1$clay + 0.022*OM - 
              0.018*(cover1$sand*OM) - 0.027*(cover1$clay*OM) -
              0.584*(cover1$sand*cover1$clay) + 0.078

theta_S_33 = theta_S_33t + (0.636*theta_S_33t - 0.107)

theta_S = theta_33 + theta_S_33 - 0.097*cover1$sand + 0.043

# Add to cover
cover1$WHC <- theta_33

################################################################################
# Aside: Calculate WHC for very coarse and very fine soils

sandx <- c(0.9,0.1)
clayx <- c(0.05, 0.75)

theta_33t = -0.251*sandx + 0.195*clayx + 0.011*OM + 
  0.006*(sandx*OM) - 0.027*(clayx*OM) +
  0.452*(sandx*clayx) + 0.299

theta_33 = theta_33t + (1.283*(theta_33t^2) - 0.374*theta_33t - 0.015)

theta_S_33t = 0.278*sandx + 0.034*clayx + 0.022*OM - 
  0.018*(sandx*OM) - 0.027*(clayx*OM) -
  0.584*(sandx*clayx) + 0.078

theta_S_33 = theta_S_33t + (0.636*theta_S_33t - 0.107)

theta_S = theta_33 + theta_S_33 - 0.097*sandx + 0.043

theta_33

################################################################################

# Set up colors to display dominance
cols <- rev(c("#40004bb3","#40004bb3","#762a83b3","#9970abb3","#c2a5cfb3","#808080B3",
              "#a6dba0b3","#5aae61b3","#1b7837b3","#00441bb3"))
# Make cuts
cover1$dom_bin <- cut(cover1$dominance, breaks = c(seq(-1,0,by = 0.25),0.0000001,seq(0,1,by=0.25)[-1],1.1))


# Plot with WHC to mimic Sala et al. 1997
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(cover1$CORR_PT~c(cover1$WHC*-1), pch = 16, 
     col = cols[cover1$dom_bin],
     cex = 2, #cex = 0.01*plots$map, 
     xlim = c(-0.44,-0.12), ylim = c(-1,1))


cor.test(cover1$dominance, cover1$CORR_PT)
# r = -0.36, p = 0.009
cor.test(cover1$dominance, cover1$WHC)
# r = -0.08, p = 0.595
cor.test(cover1$dominance, cover1$sand)
# r = 0.13, p = 0.351
cor.test(cover1$dominance, cover1$clay)
# r = -0.024, p = 0.870
cor.test(cover1$dominance, plots$map)
# r = -0.16, p = 0.266
cor.test(cover1$dominance, plots$mat)
# r = 0.422, p = 0.002

summary(lm(cover1$dominance~cover1$WHC+cover1$CORR_PT))
# R2 = 0.13, p = 0.033

summary(lm(cover1$dominance~cover1$WHC*cover1$CORR_PT))
# R2 = 0.14, p = 0.073

summary(lm(cover1$dominance~cover1$sand+cover1$CORR_PT))

table(cover1$plot == plots$plot)

summary(lm(cover1$dominance~cover1$CORR_PT+plots$map))


# Save to file
png(file.path(figdir,"Chapter1_Fig4a.png"), width = 4, height = 4, units = "in", res = 300)
par(mar= c(3,3,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(cover1$CORR_PT~c(cover1$WHC*-1), pch = 16, 
     col = cols[cover1$dom_bin],
     cex = 2,#0.01*plots$map,
     xlim = c(-0.55, -0.1), ylim = c(-1,1),
     ylab =  "", 
     xlab = "", xaxt = "n", yaxt = "n")
mtext("Soil texture", 1, line = 2, cex = 1.5)
mtext("PPT-Temp overlap", 2, line = 1.5, cex = 1.5)
axis(1, at = c(-0.52,-0.12), labels = c("fine","coarse"), tick = FALSE,
     cex.axis = 1.2)
axis(2, at = c(-0.85, 0.85), labels = c("low","high"), tick = FALSE,
     cex.axis = 1.2)
#text(expression(R^2 ~ "= 0.13"), x = 0.5, y = 0.97)
#text("p = 0.032", x = 0.5, y = 0.84)
dev.off()


# With Sand
plot(cover1$CORR_PT~cover1$sand, pch = 16, 
     col = cols[cover1$dom_bin],
     cex = 2,#0.005*plots$map, 
     xlim = c(0.15,0.85), ylim = c(-1,1))

# Use stepwise regression with WHC, CORR_PT, and map to predict dominance
mnull <- lm(cover1$dominance ~ 1)
mfull <- lm(cover1$dominance ~ cover1$WHC*cover1$CORR_PT*plots$map)
mclim <- step(mnull, scope=list(upper=mfull),direction = "both",
                 trace = 1)
summary(mclim)
par(mfrow = c(2,2), mgp = c(1,0.1,0), tcl = 0.1, mar = c(2,2,2,1))
plot(mclim)
# R2 = 0.1904, p = 0.006

# Plot dominance
cols_soil <- cut(cover1$WHC, breaks = c(0.38,0.4,0.42,0.44,0.46,0.48,0.5))
par(mfrow=c(1,1))
plot(cover1$dominance~cover1$CORR_PT, pch = 16, 
     col = viridis::viridis(6, alpha = 0.5)[cols_soil],
     cex = 0.005*plots$map, 
     xlim = c(-1,1), ylim = c(-1,1))
abline(h=0, lty = 2)
abline(v=0, lty = 2)

cols_sand <- cut(cover1$sand, breaks = seq(0.2,0.9,by=0.1))
par(mfrow=c(1,1))
plot(cover1$dominance~cover1$CORR_PT, pch = 16, 
     col = viridis::viridis(7, alpha = 0.5)[cols_sand],
     cex = 2,#0.005*plots$map, 
     xlim = c(-1,1), ylim = c(-1,1))
abline(h=0, lty = 2)
abline(v=0, lty = 2)


cor.test(cover1$dominance, cover1$CORR_PT)
# r = -0.3607, p = 0.0093
cor.test(cover1$dominance, cover1$WHC)
# r = -0.08, p = 0.595
cor.test(cover1$dominance, cover1$sand)
# r = 0.1333, p = 0.3507

# Make a map of dominance across the West
# Load sagebrush region
art <- rast('C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/01_Data/CoreARTR_combined_DayMet_cropped&trimmed.tif')

# Make file
png(file.path(figdir,"Sites_in_SageRegion_dominance.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(1,1,1,1), tcl = -0.2)
maps::map("state", regions = c("Wyoming","Nevada","Idaho","Colorado","Arizona",
                               "Utah","Washington","Oregon","Montana","California",
                               "New Mexico"),
          xlim = c(-125,-101), ylim = c(30,50))
image(art, col = rgb(0.5,0.5,0.5,0.5),add = T)
points(plots$latitude~plots$longitude, pch = 16, cex=1.5,,
       col = cols[cover1$dom_bin])
points(plots$latitude~plots$longitude, pch = 1, cex=1.5)
maps::map.scale(x = -124, y = 32, relwidth = 0.25, ratio = FALSE, cex = 0.7)
GISTools::north.arrow(x = -102, y = 46.5, lab = "N", len = 0.6, col = "black")
dev.off()

png(file.path(figdir,"Sites_in_SageRegion_dominance_and_MAP_legend.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(1,1,1,1), las = 0, mgp = c(1,0.1,0))
plot(c(1:10)~1, pch = 16, col = "white", yaxt = "n", xaxt = "n", bty = "n",
     xlab = "", ylab ="")
for (i in 1:(length(cols)-1)){
  polygon(c(5,5,5.5,5.5), c(i,i+1,i+1,i), col = cols[i], lwd = 2)
}
axis(side = 2, at = c(2,9), labels = c("Grass","Shrub"), pos = 5, tick = F,
     cex.axis = 1.7)
dev.off()

# Show points in climate space
png(file.path(figdir,"ClimateSpace&dominance.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(2,2,1,1), tcl = 0.1, mgp = c(1,0.1,0))
plot(plots$map~plots$mat, pch = 16, cex=2,
       col = cols[cover1$dom_bin],
     xlab = "Mean annual temperature (C)",
     ylab = "Mean annual precipitation (mm)")
points(plots$map~plots$mat, pch = 1, cex=2, lwd = 2)
dev.off()



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

# Visualize 
table(cover1$plot == results$plot)
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(results$avg_trans_0_200cm~results$prop_shallow, pch = 16, 
     col = cols[cover1$dom_bin],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(35, 300))


# Use stepwise regression with prop_shallow & totaltransp to predict dominance
mnull <- lm(cover1$dominance ~ 1)
mfull <- lm(cover1$dominance ~ results$prop_shallow*results$avg_trans_0_200cm)
msw <- step(mnull, scope=list(upper=mfull),direction = "both",
              trace = 1)
summary(msw)
par(mfrow = c(2,2), mgp = c(1,0.1,0), tcl = 0.1, mar = c(2,2,2,1))
plot(msw)
# R2 = 0.2974, p = 0.0002

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

#Reorder
dung2 <- dung2a[order(dung2a$plot),]


# Visualize 
table(cover1$plot == dung2$plot)
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(results$avg_trans_0_200cm~results$prop_shallow, pch = 16, 
     col = cols[cover1$dom_bin],
     cex = 20*(dung2$cow+0.05),
     xlim = c(0.3,0.7), ylim = c(35, 300))

# Use stepwise regression with prop_shallow & totaltransp & dung to predict dominance
mnull <- lm(cover1$dominance ~ 1)
mfull <- lm(cover1$dominance ~ results$prop_shallow*results$avg_trans_0_200cm*dung2$cow*dung2$horse)
msw_cow<- step(mnull, scope=list(upper=mfull),direction = "both",
              trace = 1)
summary(msw_cow)
par(mfrow = c(2,2), mgp = c(1,0.1,0), tcl = 0.1, mar = c(2,2,2,1))
plot(msw_cow)
# R2 = 0.3743, p <0.001

# Use stepwise regression with climate, sand, & dung to predict dominance
mnull <- lm(cover1$dominance ~ 1)
mfull <- lm(cover1$dominance ~ cover1$CORR_PT*plots$map*dung2$cow*dung2$horse*cover1$sand)
mclim_cow<- step(mnull, scope=list(upper=mfull),direction = "both",
               trace = 1)
summary(mclim_cow)
par(mfrow = c(2,2), mgp = c(1,0.1,0), tcl = 0.1, mar = c(2,2,2,1))
plot(mclim_cow)
# R2 = 0.26, p = 0.0007

# Check correlation between dominance/PG and cow and horse
cor.test(cover1$dominance, dung2$cow)
# r = -0.37, p = 0.007
cor.test(cover1$dominance, dung2$horse)
# r = 0.10, p = 0.502
cor.test(cover1$PG, dung2$cow)
# r = 0.3745, p = 0.0068
cor.test(cover1$PG, dung2$horse)
# r = -0.18, p = 0.214
cor.test(cover1$SH, dung2$cow)
# r = 0.01, p = 0.964
cor.test(cover1$SH, dung2$horse)
# r = -0.22, p = 0.116
cor.test(results$avg_trans_0_200cm, dung2$cow)
# r = 0.3130, p = 0.0253
cor.test(results$avg_trans_0_200cm, dung2$horse)
# r = -0.27, p = 0.0553
cor.test(plots$map, dung2$cow)
# r = 0.22, p = 0.1174
cor.test(plots$map, dung2$horse)
# r = -0.22, p = 0.122
cor.test(dung2$cow, dung2$horse)
# r = -0.2, p = 0.159

cor.test(plots$mat, dung2$cow)
# r = -0.31, p = 0.028
cor.test(plots$mat, dung2$horse)
# r = 0.14, p = 0.339


################################################################################
# Step 5: Create figures relating climate and soilwater to dominance

# A) MAP

png(file.path(figdir, "Dominance_vs_map.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, plots$map)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~plots$map, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "Mean annual precipitation (mm)", cex = 1.4, ylim = c(-1,1))
points(cover1$dominance~plots$map, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# B)  Total transpiration

png(file.path(figdir, "Dominance_vs_totaltranspiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, results$avg_trans_0_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~results$avg_trans_0_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "Mean annual transpiration (mm)", cex = 1.4, ylim = c(-1,1))
points(cover1$dominance~results$avg_trans_0_200cm, lwd = 1.7, col = 'black', cex = 1.4)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# C)  Shallow proportion of transpiration

png(file.path(figdir, "Dominance_vs_proportion_shallow_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, results$prop_shallow)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~results$prop_shallow, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "",
     cex = results$avg_trans_0_200cm*0.015, 
     xlim = c(0.3,0.7),ylim = c(-1,1))
mtext("Proportion of total transpiration from\nshallow soils (0-30 cm)", 1, line = 2)
points(cover1$dominance~results$prop_shallow, lwd = 1.7, col = 'black', cex = results$avg_trans_0_200cm*0.015)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# D) Total shallow soil transpiration

png(file.path(figdir, "Dominance_vs_total_shallow_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, results$avg_trans_0_30cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~results$avg_trans_0_30cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "Mean shallow soil transpiration (0-30 cm)", cex = 2, ylim = c(-1,1))
points(cover1$dominance~results$avg_trans_0_30cm, lwd = 1.7, col = 'black', cex = 2)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# E) Total deep soil transpiration

png(file.path(figdir, "Dominance_vs_total_deep_transpiration.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, results$avg_trans_30_200cm)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~results$avg_trans_30_200cm, pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "Mean deep soil transpiration (30-200 cm)", cex = 2, ylim = c(-1,1))
points(cover1$dominance~results$avg_trans_30_200cm, lwd = 1.7, col = 'black', cex = 2)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

# F) Shrub cover

png(file.path(figdir, "Dominance_vs_shrub_cover.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, cover1$SH)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~c(cover1$SH*100), pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "Shrub cover (%)", cex = results$avg_trans_0_200cm*0.015, ylim = c(-1,1))
points(cover1$dominance~c(cover1$SH*100), lwd = 1.7, col = 'black', cex = results$avg_trans_0_200cm*0.015)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p = ",round(pg_cor$p.value, 3)),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

png(file.path(figdir, "Dominance_vs_Pgrass_cover.png"), width = 5, height = 4.5, units = "in",
    res = 300)
pg_cor <- cor.test(cover1$dominance, cover1$PG)
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(cover1$dominance~c(cover1$PG*100), pch = 16, col = rgb(0.5,0.5,0.5,0.7),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "Perennial grass cover (%)", cex = results$avg_trans_0_200cm*0.015, ylim = c(-1,1))
points(cover1$dominance~c(cover1$PG*100), lwd = 1.7, col = 'black', cex = results$avg_trans_0_200cm*0.015)
legend("bottomright", legend = paste0("r = ",round(pg_cor$estimate,2),
                                      ", p < 2.2e-16"),
       bty = "n")
abline(h = 0, lty = 2)
dev.off()

################################################################################
# Step 6: Final figure showing dominance ~ amount + depth distribution of soil water

# Create binary variable for wet vs dry sites
mean(results$avg_trans_0_200cm)
results$dry <- as.numeric(results$avg_trans_0_200cm < mean(results$avg_trans_0_200cm))

# Create a second binary variable for wet vs. dry sites based on median
results$dry_median <- as.numeric(results$avg_trans_0_200cm < median(results$avg_trans_0_200cm))

# Now create variable with several levels (1 = very dry, 2 = moderate, 3 = very wet)
results$moisture <- NA
results[results$avg_trans_0_200cm < quantile(results$avg_trans_0_200cm, 0.33),]$moisture <- 1
results[results$avg_trans_0_200cm >= quantile(results$avg_trans_0_200cm, 0.66),]$moisture <- 2
results[results$avg_trans_0_200cm >= quantile(results$avg_trans_0_200cm, 0.33) &
          results$avg_trans_0_200cm < quantile(results$avg_trans_0_200cm, 0.66)  ,]$moisture <- 3
# Set as a factor
results$moisture <- as.factor(results$moisture)

summary.aov(lm(cover1$dominance~results$prop_shallow+results$moisture))
summary(lm(cover1$dominance~results$prop_shallow+results$moisture))

summary(lm(cover1$dominance~results$prop_shallow+results$dry))

summary(lm(cover1$dominance~results$prop_shallow+results$dry_median))

summary(lm(cover1$dominance~results$prop_shallow+results$avg_trans_0_200cm))


# Look at 3 moisure levels
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(cover1$dominance~results$prop_shallow, pch = 16, 
     col = c(rgb(1,0.1,0,0.5), rgb(0.5,0.5,0.5,0.5), rgb(0,0.1,1,0.5))[results$moisture],
     cex = 20*(dung2$cow+0.05),
     xlim = c(0.3,0.7), ylim = c(-1,1))
abline(h=0, lty = 2)
moist1 <- smooth.spline(x = results[results$moisture == 1, ]$prop_shallow, 
                        y = cover1[results$moisture == 1, ]$dominance,
                        all.knots = FALSE, nknots = 10)
lines(moist1, col = "darkred", lwd=2)
moist2 <- smooth.spline(x = results[results$moisture == 2, ]$prop_shallow, 
                        y = cover1[results$moisture == 2, ]$dominance,
                        all.knots = FALSE, nknots = 10)
lines(moist2, col = "darkgrey", lwd=2)
moist3 <- smooth.spline(x = results[results$moisture == 3, ]$prop_shallow, 
                        y = cover1[results$moisture == 3, ]$dominance,
                        all.knots = FALSE, nknots = 10)
lines(moist3, col = "darkblue", lwd=2)


# Look at just wet and dry plots
par(mar= c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(cover1$dominance~results$prop_shallow, pch = 16, 
     col = c(rgb(0,0.1,1,0.5), rgb(1,0.1,0,0.5))[results$dry+1],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(-1,1))
abline(h=0, lty = 2)
dry1 <- smooth.spline(x = results[results$dry == 1, ]$prop_shallow, 
                      y = cover1[results$dry == 1, ]$dominance,
                      all.knots = FALSE, nknots = 10, spar = 0.7)
lines(dry1, col = "darkred", lwd=2)
dry2 <- smooth.spline(x = results[results$dry == 0, ]$prop_shallow, 
                      y = cover1[results$dry == 0, ]$dominance,
                      all.knots = FALSE, nknots = 10, spar = 0.7)
lines(dry2, col = "darkblue", lwd=2)

# Now with linear models only
m <- lm(cover1$dominance~results$prop_shallow+results$dry)

png(file.path(figdir,"Chapter1_Fig4b.png"), width = 4, height = 4, units = "in", res = 300)
par(mar= c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(cover1$dominance~results$prop_shallow, pch = 16, 
     col = c("#01665eb3", "#8c510ab3")[results$dry+1],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(-1,1),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "")
#points(cover1$dominance~results$prop_shallow, cex = 2, lwd = 2)
mtext("Proportion of total transpiration from\nshallow soils (0-30 cm)", 1, line = 2)
abline(h=0, lty = 2)
curve(m$coefficients[1] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#01665e") # wet plots
curve(m$coefficients[1] + m$coefficients[3] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#8c510a") # dry plots
mtext("             Dry: Y = 1.65 - 2.56*X", side = 1, line = -2.2, adj = 0)
mtext("             Wet: Y = 1.39 - 2.56*X", side = 1, line = -1, adj = 0)
# Wet
points(x = 0.325, y = -1.009, cex = 2, col = "#01665eb3", pch = 16)
#points(x = 0.325, y = -1.009, cex = 2, lwd = 2)
lines(x = c(0.3,0.35), y = rep(-1.009,2), lwd = 3, col = "#01665e")
# Dry
points(x = 0.325, y = -0.86, cex = 2, col = "#8c510ab3", pch = 16)
#points(x = 0.325, y = -0.86, cex = 2, lwd = 2)
lines(x = c(0.3,0.35), y = rep(-0.86,2), lwd = 3, col = "#8c510a")
text(expression(R^2 ~ "= 0.29"), x = 0.67, y = 0.97)
text("p < 0.001", x = 0.67, y = 0.84)
dev.off()
