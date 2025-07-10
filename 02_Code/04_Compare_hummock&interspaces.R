################################################################################
# Purpose: Evaluate hummock effects on transpiration and interactions
#
#
# Rachel R. Renne
# April 6, 2025
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
plot(cover1$PG ~ cover1$SH, pch = 16, col = rgb(0.5,0.5,0.5,0.7), cex = 2)

# Calculate shrub/grass dominance index:
cover1$dominance <- (cover1$SH - cover1$PG)/(cover1$SH + cover1$PG)


################################################################################
# Step 2: Get RIIs for PG

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

# Pull out individual ftypes for comparison
pgrass <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype == "PG",]

################################################################################
# Step 3: Get Water Balance variables

# Set up project directory
dir_int <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/04_Simulation_Experiment/Interspace_simulations/Interspace_simulations_20250114"
dir_hum <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/04_Simulation_Experiment/Hummock_simulations/Hummock_simulations_20250114"

# Create a plotlist:
plotlist1 <- read.csv(file.path(dir_int,"1_Input","SWRuns_InputMain_interspace_simulations_v12.csv"))
plotlist <- plotlist1$Label

# Set up for results
results <- data.frame(plot = plotlist, ievap_surf = NA, ievap_soil = NA,
                      itrans_shallow = NA, itrans_deep = NA, idrainage = NA,
                      hevap_surf = NA, hevap_soil = NA,
                      htrans_shallow = NA, htrans_deep = NA, hdrainage = NA)

# Database directory
dir_dat_int <- file.path(dir_int, "4_Simulation")
dir_dat_hum <- file.path(dir_hum, "4_Simulation")

for (dir_dat in c(dir_dat_int, dir_dat_hum)){
  
  # Set up dir_prj
  if (dir_dat == dir_dat_int){
    dir_prj <- dir_int
  } else { dir_prj <- dir_hum}
  
  # Database name
  outDB_fname <- file.path(dir_dat, "dbOutput.sqlite3")
  
  # Connect to output database
  outDB_con <- RSQLite::dbConnect(RSQLite::SQLite(), outDB_fname)
  
  # Get data for entire tables (daily aggregations)
  
  # Average daily transpiration
  trans <- get.Table_Scenario(outDB_fname,
                              responseName = "aggregation_doy_Transpiration", MeanOrSD = "Mean",
                              scenario = "Current", header = TRUE)
  # Average deep drainage
  deep <- get.Table_Scenario(outDB_fname,
                             responseName = "aggregation_doy_DeepDrainage", MeanOrSD = "Mean",
                             scenario = "Current", header = TRUE)
  
  # Evaporation from soil (bare soil evaporation)
  esoil <- get.Table_Scenario(outDB_fname,
                              responseName = "aggregation_doy_EvaporationSoil", MeanOrSD = "Mean",
                              scenario = "Current", header = TRUE)
  
  # Evaporation from surface (evaporation of intercepted rain)
  esurface <- get.Table_Scenario(outDB_fname,
                                 responseName = "aggregation_doy_EvaporationSurface", MeanOrSD = "Mean",
                                 scenario = "Current", header = TRUE)
  
  # Clean up
  RSQLite::dbDisconnect(outDB_con)  
  
  ##############################################################################
  # Populate results
  
  for (i in 1:length(plotlist)){
    
    # Get surface evaporation
    if (dir_prj == dir_int){
      # Surface evaporation
      results$ievap_surf[i] <- sum(esurface[esurface$Label == plotlist[i],51:416])
      
      # Soil evaporation
      results$ievap_soil[i] <- sum(esoil[esoil$Label == plotlist[i],51:416])
      
      # Shallow transpiration
      results$itrans_shallow[i] <- sum(trans[trans$Label == plotlist[i],52:417][1:4,])
      
      # Deep transpiration
      results$itrans_deep[i] <- sum(trans[trans$Label == plotlist[i],52:417][5:15,])
      
      # Deep drainage
      results$idrainage[i] <- sum(deep[deep$Label == plotlist[i],51:416])
      
    } else {
      # Surface evaporation
      results$hevap_surf[i] <- sum(esurface[esurface$Label == plotlist[i],51:416])
      
      # Soil evaporation
      results$hevap_soil[i] <- sum(esoil[esoil$Label == plotlist[i],51:416])
      
      # Shallow transpiration
      results$htrans_shallow[i] <- sum(trans[trans$Label == plotlist[i],52:417][1:4,])
      
      # Deep transpiration
      results$htrans_deep[i] <- sum(trans[trans$Label == plotlist[i],52:417][5:15,])
      
      # Deep drainage
      results$hdrainage[i] <- sum(deep[deep$Label == plotlist[i],51:416])
    }
  }
}
    
# Separate results between hummock and interspace
ires <- results[,c(1:6)]
names(ires)[2:6] <- c("evap_surf","evap_soil","trans_shallow","trans_deep","drainage")
hres <- results[,c(1,7:11)]
names(hres)[2:6] <- c("evap_surf","evap_soil","trans_shallow","trans_deep","drainage")

################################################################################
# Step 4: Get hummock data

# Read in hummock data
hum <- read.csv(file.path(datadir,"hummocks.csv"))

# Merge ghum and hum
hum1 <- merge(hum, plotlist1[,1:2], by.x = "plot", by.y = "Label", all = TRUE)

# Check for missing data
hum1[is.na(hum1$site_id),]

# Get average depth across all 5 samples
hum3 <- data.frame(plot = sort(unique(hum1$plot)),
                   depth = tapply(hum1$mineralsoil, hum1$plot, function(x){mean(x, na.rm = TRUE)}))

# Now round depths (and set negative depth to 0)
hum3[hum3$depth < 0,]$depth <- 0
hum3$depth <- round(hum3$depth)

# Read in soil texture
texture <- read.csv(file.path(datadir,"soiltexture.csv"))

# Get interspace surface texture
itexture <- texture[!is.na(texture$depth_0) & texture$depth_0 == 0,c(1,15,16)]

table(hum3$plot == itexture$plot)

# Add onto hum3
hum3$isand <- itexture$sand
hum3$iclay <- itexture$clay

# Get hummock average texture
htexture <- data.frame(hsand = tapply(texture[!is.na(texture$transect),"sand"],texture[!is.na(texture$transect),"plot"], mean),
                       hclay = tapply(texture[!is.na(texture$transect),"clay"],texture[!is.na(texture$transect),"plot"], mean))
htexture$plot <- row.names(htexture)         

# Add onto hum3
hum4 <- merge(hum3, htexture, by = "plot", all = T)

table(hres$plot == hum4$plot)
table(ires$plot == hum4$plot)
table(pgrass$plot == hum4$plot)
table(plots$plot == hum4$plot)
table(hum4$plot == cover1$plot)
# All true

################################################################################
# Step 5: Relate hummock and interspace characteristics to interactions

hres$total_trans <- apply(hres[,5:6],1,sum)
ires$total_trans <- apply(ires[,5:6],1,sum)

hres$prop_shallow <- hres$trans_shallow/(hres$trans_shallow+hres$trans_deep)
ires$prop_shallow <- ires$trans_shallow/(ires$trans_shallow+ires$trans_deep)


# Interaction is ONLY related to sand content, but not difference in sand content
# Probably has to do with how sand favors shrubs a bit more
cor.test(pgrass$pcomm_cover_rii, hum4$depth)
# r = 0.169, p = 0.245
cor.test(pgrass$pcomm_cover_rii, c(hum4$hsand - hum4$isand))
# r = 0.074, p = 0.618
cor.test(pgrass$pcomm_cover_rii, hum4$hsand)
# r = 0.363, p = 0.011
cor.test(pgrass$pcomm_cover_rii, hum4$isand)
# r = 0.302, p = 0.035
cor.test(pgrass$pcomm_cover_rii, c(hres$total_trans - ires$total_trans))
# r = -0.01, p = 0.964
cor.test(pgrass$pcomm_cover_rii, c(hres$evap_surf - ires$evap_surf))
# r = -0.05, p = 0.735
cor.test(pgrass$pcomm_cover_rii, c(hres$evap_soil - ires$evap_soil))
# r = 0.012, p = 0.934
cor.test(pgrass$pcomm_cover_rii, c(hres$drainage - ires$drainage))
# r = -0.038, p = 0.798

# Look at relationships between hummocks and water
cor.test(hum4$depth, c(hum4$hsand - hum4$isand))
# r = 0.349, p = 0.013 # Deeper hummocks have larger sand diff
cor.test(hum4$depth, c(hres$total_trans - ires$total_trans))
# r = 0.20, p = 0.16 -- depth doesn't matter
cor.test(c(hum4$hsand-hum4$isand), c(hres$total_trans - ires$total_trans))
# r = 0.329, p = 0.02 # But sand difference matters...
cor.test(hum4$hsand, c(hres$total_trans - ires$total_trans))
# r = 0.156, p = 0.281 # Sand content doesn't matter
cor.test(hum4$isand, c(hres$total_trans - ires$total_trans))
# r = -0.026, p = 0.854]
cor.test(hum4$depth, ires$total_trans)
# r = -0.223, p = 0.116
cor.test(hum4$depth, hres$total_trans)
# r = -0.211, p = 0.137
cor.test(hum4$depth, plots$map)
# r = -0.34, p = 0.015 # Wetter sites = smaller hummocks
cor.test(hum4$depth, plots$mat)
# r = 0.289, p = 0.04 # Hotter sites = taller hummocks
cor.test(c(hum4$hsand-hum4$isand), plots$map)
# r = -0.352, p = 0.012 # Wetter sites = smaller sand difference
cor.test(c(hum4$hsand-hum4$isand), plots$mat)
# r = 0.26, p = 0.063 # Hotter sites = larger sand difference

t.test(x = hum4$hsand, y = hum4$isand, paired = TRUE)

# But MAP and MAT aren't related to diff's in transpiration between hum & int
cor.test(plots$mat, c(hres$total_trans - ires$total_trans))
# r = 0.056, p = 0.694
cor.test(plots$map, c(hres$total_trans - ires$total_trans))
# r = 0.041, p = 0.773
cor.test(c(hres$prop_shallow-ires$prop_shallow), c(hres$total_trans - ires$total_trans))
# r = -0.466, p < 0.001 # Change in prop_shallow = increase 

# How about dominance--not related to ANY of this stuff
cor.test(cover1$dominance, hum4$depth)
# r = 0.034, p = 0.815
cor.test(cover1$dominance, c(hum4$hsand - hum4$isand))
# r = 0.137, p = 0.344
cor.test(cover1$dominance, hum4$hsand)
# r = 0.221, p = 0.123
cor.test(cover1$dominance, hum4$isand)
# r = 0.137, p = 0.338
cor.test(cover1$dominance, c(hres$total_trans - ires$total_trans))
# r = 0.182 p = 0.200
cor.test(cover1$dominance, c(hres$evap_surf - ires$evap_surf))
# r = -0.199, p = 0.163
cor.test(cover1$dominance, c(hres$evap_soil - ires$evap_soil))
# r = -0.20, p = 0.163
cor.test(cover1$dominance, c(hres$drainage - ires$drainage))
# r = 0.106, p = 0.460


plot(pgrass$pcomm_cover_rii ~ c(hres$total_trans - ires$total_trans),
     cex = 0.5*(hum4$depth+1), pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     ylim = c(-1,1))
abline(v = 0, lty = 2)

plot(pgrass[hum4$depth > 0, ]$pcomm_cover_rii ~ c(hum4[hum4$depth > 0, ]$hsand - hum4[hum4$depth > 0, ]$isand),
     cex = 0.5*(hum4[hum4$depth > 0, ]$depth), pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     ylim = c(-1,1))
abline(v = 0, lty = 2)

plot(c(hres$total_trans - ires$total_trans)~hum4$depth,
     pch = 16, col = rgb(0.5,0.5,0.5,0.5), cex = 2)
abline(h=0, lty = 2)

plot(pgrass$pcomm_cover_rii ~ c(hres$total_trans - ires$total_trans),
     pch = 16, col = rgb(0.5,0.5,0.5,0.5), ylim = c(-1,1),
     cex = 0.5*(hum4$depth+1))
abline(h=0, lty = 2)

par(mar = c(2,2,1,1))
plot(c(hres$total_trans - ires$total_trans)~ires$prop_shallow,
     pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     cex = 0.5*(hum4$depth))
abline(h=0, lty = 2)

plot(c(hres$total_trans - ires$total_trans)~plots$mat,
     pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     cex = 0.5*(hum4$depth+0.5))
abline(h=0, lty = 2)

t.test(x = hres$total_trans, y = ires$total_trans, paired = TRUE)
t.test(x = hres$prop_shallow, y = ires$prop_shallow, paired = TRUE)

summary(hres$total_trans-ires$total_trans)
summary(hres$prop_shallow-ires$prop_shallow)

summary(lm(c(hres$total_trans - ires$total_trans)~c(hum4$hsand*100-hum4$isand*100)))


summary(lm(c(hres$total_trans - ires$total_trans)~c(hum4$hsand*100-hum4$isand*100)+c(hres$prop_shallow*100 - ires$prop_shallow*100)))

plot(c(hres$total_trans - ires$total_trans)~c(hum4$hsand*100-hum4$isand*100),
     pch = 16, col = rgb(0.5,0.5,0.5,0.5), cex = plots$map*0.01)
abline(h=0,lty=2)

plot(c(hres$total_trans - ires$total_trans)~c(hres$prop_shallow*100 - ires$prop_shallow*100),
     pch = 16, col = rgb(0.5,0.5,0.5,0.5), cex = plots$map*0.01)
abline(h=0,lty=2)

plot(plots$map~plots$mat,
     pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     cex = 0.5*(hum4$depth+1))

################################################################################
# Step 6: Make water balance comparison figure

# Relativize the water balance
rel_hres <- apply(hres[,2:6],1,function(x){x/sum(x)})
rel_ires <- apply(ires[,2:6],1,function(x){x/sum(x)})

# Calculate average of each
avg_hres <- apply(rel_hres[,hum4$depth > 0],1, mean)
avg_ires <- apply(rel_ires[,hum4$depth > 0],1, mean)

# Make stacked bar plot of each
par(mar = c(2,1,1,1), mgp = c(1,0.1,0), tcl = 0.1)
barplot(cbind(rev(avg_hres), rev(avg_ires)),
        col = c("maroon","darkblue","blue","darkorange","darkred"),
        names = c("Hummock","Interspace"), ylim = c(0,1),
        border = NA, space = 0.01, yaxt = "n")

# Try histograms of differences instead
png(file.path(figdir, "WaterBalance_hummock_interspace.png"), width = 4, height = 5, units = "in",
    res = 300)
par(mgp = c(1,0.1,0), tcl = 0.1)
layout.matrix <- matrix(c(1:6), nrow = 6, byrow = TRUE)
layout(layout.matrix, widths = rep(1,6), heights = c(rep(1,5),0.5))
layout.show(6)

par(mar = c(0,2,0.25,1))
hist(hres$evap_surf - ires$evap_surf, xlim = c(-17.5,17.5),
     breaks = seq(-17.5,17.5, by = 2.5),
     main = "", xaxt="n", ylab = "Freq", ylim = c(0,25),
     cex.axis = 1.2, cex.lab = 1.2)
mtext("Surface evaporation  ", side = 3, line = -1.25, adj = 1)
box()
abline(h=0)
abline(v = 0, lty = 2, lwd = 1.5)

hist(hres$evap_soil - ires$evap_soil, xlim = c(-17.5,17.5),
     breaks = seq(-17.5,17.5, by = 2.5),
     main = "", xaxt="n", ylab = "Freq", ylim = c(0,25),
     cex.axis = 1.2, cex.lab = 1.2)
mtext("Bare soil evaporation  ", side = 3, line = -1.25, adj = 1)
box()
abline(h=0)
abline(v = 0, lty = 2, lwd = 1.5)

hist(hres$trans_shallow - ires$trans_shallow, xlim = c(-17.5,17.5),
     breaks = seq(-17.5,17.5, by = 2.5),
     main = "", xaxt="n", ylab = "Freq", ylim = c(0,25),
     cex.axis = 1.2, cex.lab = 1.2)
mtext("Shallow transpiration  ", side = 3, line = -1.25, adj = 1)
box()
abline(h=0)
abline(v = 0, lty = 2, lwd = 1.5)

hist(hres$trans_deep - ires$trans_deep, xlim = c(-17.5,17.5),
     breaks = seq(-17.5,17.5, by = 2.5),
     main = "", xaxt="n", ylab = "Freq", ylim = c(0,25),
     cex.axis = 1.2, cex.lab = 1.2)
mtext("Deep transpiration  ", side = 3, line = -1.25, adj = 1)
box()
abline(h=0)
abline(v = 0, lty = 2, lwd = 1.5)

hist(hres$drainage - ires$drainage, xlim = c(-17.5,17.5),
     breaks = seq(-17.5,17.5, by = 2.5),
     main = "", xaxt="n", ylab = "Freq", ylim = c(0,25),
     cex.axis = 1.2, cex.lab = 1.2)
mtext("Drainage  ", side = 3, line = -1.25, adj = 1)
box()
abline(h=0)
abline(v = 0, lty = 2, lwd = 1.5)

par(mgp = c(-1,-1.5,-0), tcl = 0.1, mar= c(0,2,0,1))
plot(1, xlim = c(-17.5,17.5), xaxt = "n",yaxt = "n", 
     col = "white", ylab = "", xlab = "", bty = "n")
axis(side = 3, at = seq(-15,15,by = 5), labels = seq(-15,15,by=5), cex.axis = 1.5)
mtext("Hummock - Interspace (mm)", side = 1, line = -1, cex = 1.2)
dev.off()

par(mfrow = c(1,1), mgp = c(1,0.1,0), mar = c(2,2,1,1))
cor_hum <- cor.test(hum4$depth, c(hum4$hsand-hum4$isand))
# r = 0.349, p = 0.013
plot(hum4$depth~c(hum4$hsand-hum4$isand))

