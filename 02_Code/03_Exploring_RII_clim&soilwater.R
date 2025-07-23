################################################################################
# Purpose: Evaluate shrub-grass interactions across 51 shrub-grass interaction 
#          sites in sagebrush habitat across the western U.S.
#
#
# Rachel R. Renne
# April 2, 2025
# Updated: July 10, 2025
################################################################################

# Load libraries
library(dplyr)

# Set up directories
datadir <- '01_Data'
codedir <- '02_Code'
figdir <- '03_Figures'

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
codes1 <- codes[,c("code","GrowthForm")]
# Change GrowthForm for NO GAPS
codes1[is.na(codes1$GrowthForm),]$GrowthForm <- "NO GAPS"
names(codes1) <- c("species","ftype")
ftype = codes1

# Calculate cover of each species
cover <- lpicover_ftypecorrected(lpi, ftype, n_layers = 5,
                                 soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST","R"))

# Make dataframe wider
cover1 <- tidyr::pivot_wider(cover, id_cols = plot, names_from = species, 
                             values_from = cover, values_fill = 0)

# Calculate shrub/grass dominance index:
cover1$dominance <- (cover1$SH - cover1$PG)/(cover1$SH + cover1$PG)

################################################################################
# Step 2: Get RIIs for PGrass

# Read in community data
pcomm <- read.csv(file.path(datadir, "community.csv"))

# Calculate both microsite and canopy condition values
ms_and_ca_calcs <- pcomm_microsite(pcomm, ftype, can_edge = FALSE)

# Pull out just microsites RII
# In this case, "canopy" microsite is canopy > 0, "interspace" is canopy = 0
microsites1 <- ms_and_ca_calcs[[1]]

# Calculate some RIIs (cover, number, size, height)
microsites1$pcomm_cover_rii <- calculate_rii(canopy=microsites1$canopy_avg_cover, 
                                             interspace=microsites1$interspace_avg_cover)

microsites1$pcomm_number_rii <- calculate_rii(canopy=microsites1$canopy_avg_number, 
                                              interspace=microsites1$interspace_avg_number)

microsites1$pcomm_size_rii <- calculate_rii(canopy=microsites1$canopy_size,
                                            interspace=microsites1$interspace_size)

microsites1$pcomm_height_rii <- calculate_rii(canopy=microsites1$canopy_height,
                                              interspace=microsites1$interspace_height)

################################################################################
# Step 3: Calculate interactions

# Pull out just perennial grasses
pgrass <- microsites1[!is.na(microsites1$ftype) & microsites1$ftype == "PG",]

# Combine with plots dataframe
pgrass <- merge(pgrass, plots, by = "plot", all = TRUE)

################################################################################
# Step 4: Check for dung effects

# Read in dung
dung <- read.csv(file.path(datadir, "dung.csv"))

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

# Check correlations
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

################################################################################
# Step 5: Calculate correlations between RIIs and 1) Climate, 2) transpiration

#####################################
# RII Pgrass cover
summary(pgrass$pcomm_cover_rii)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.76073 -0.35800 -0.03226  0.07302  0.43767  1.00000        2 
# Note: two NA's at sites with no grass in any quadrats

cor.test(pgrass$pcomm_cover_rii, pgrass$map)
# r = -0.20, p = 0.174
cor.test(pgrass$pcomm_cover_rii, pgrass$mat)
# r = 0.31, p = 0.030
cor.test(pgrass$pcomm_cover_rii, pgrass$avg_PET)
# r = 0.20, p = 0.159

cor.test(pgrass$pcomm_cover_rii, pgrass$avg_trans_0_200cm)
# r = -0.28, p = 0.049
cor.test(pgrass$pcomm_cover_rii, pgrass$avg_trans_0_30cm)
# r = -0.44, p = 0.001
cor.test(pgrass$pcomm_cover_rii, pgrass$avg_trans_30_200cm)
# r = -0.11, p = 0.470
cor.test(pgrass$pcomm_cover_rii, pgrass$prop_shallow)
# r = -0.51, p < 0.001
 
#####################################
# RII Pgrass size
summary(pgrass$pcomm_size_rii)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.68106 -0.14715  0.02500  0.06409  0.23510  0.96522        6 

# If there are no grasses in an interspace, it should be size = 0 instead of NA
pgrass[is.na(pgrass$pcomm_size_rii) & !is.na(pgrass$canopy_size),]$interspace_size <- 0
pgrass[is.na(pgrass$pcomm_size_rii) & !is.na(pgrass$canopy_size),]$pcomm_size_rii <- 1

summary(pgrass$pcomm_size_rii)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.68106 -0.12163  0.05553  0.14049  0.37853  1.00000        2 

cor.test(pgrass$pcomm_size_rii, pgrass$map)
# r = -0.15, p = 0.289
cor.test(pgrass$pcomm_size_rii, pgrass$mat)
# r = 0.30, p = 0.038
cor.test(pgrass$pcomm_size_rii, pgrass$avg_PET)
# r = 0.23, p = 0.117

cor.test(pgrass$pcomm_size_rii, pgrass$avg_trans_0_200cm)
# r = -0.25, p = 0.085
cor.test(pgrass$pcomm_size_rii, pgrass$avg_trans_0_30cm)
# r = -0.39, p = 0.005
cor.test(pgrass$pcomm_size_rii, pgrass$avg_trans_30_200cm)
# r = -0.09, p = 0.558
cor.test(pgrass$pcomm_size_rii, pgrass$prop_shallow)
# r = -0.44, p = 0.001

#####################################
# RII Pgrass number
summary(pgrass$pcomm_number_rii)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.43243 -0.25000 -0.06897  0.02759  0.14286  1.00000        2  
# Note: two NA's at sites with no grass in any quadrats

cor.test(pgrass$pcomm_number_rii, pgrass$map)
# r = -0.28, p = 0.053
cor.test(pgrass$pcomm_number_rii, pgrass$mat)
# r = 0.34, p = 0.016
cor.test(pgrass$pcomm_number_rii, pgrass$avg_PET)
# r = 0.23, p = 0.117

cor.test(pgrass$pcomm_number_rii, pgrass$avg_trans_0_200cm)
# r = -0.30, p = 0.033
cor.test(pgrass$pcomm_number_rii, pgrass$avg_trans_0_30cm)
# r = -0.4, p = 0.005
cor.test(pgrass$pcomm_number_rii, pgrass$avg_trans_30_200cm)
# r = -0.19, p = 0.199
cor.test(pgrass$pcomm_number_rii, pgrass$prop_shallow)
# r = -0.39, p = 0.006

#####################################
# RII Pgrass height
summary(pgrass$pcomm_height_rii)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.04083  0.07765  0.14087  0.16186  0.21877  0.77887        6 

# If there are no grasses in an interspace, it should be size = 0 instead of NA
pgrass[is.na(pgrass$pcomm_height_rii) & !is.na(pgrass$canopy_height),]$interspace_height <- 0
pgrass[is.na(pgrass$pcomm_height_rii) & !is.na(pgrass$canopy_height),]$pcomm_height_rii <- 1

summary(pgrass$pcomm_height_rii)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -0.04083  0.08015  0.15853  0.23028  0.22670  1.00000        2 

cor.test(pgrass$pcomm_height_rii, pgrass$map)
# r = -0.32, p = 0.026
cor.test(pgrass$pcomm_height_rii, pgrass$mat)
# r = 0.40, p = 0.004
cor.test(pgrass$pcomm_height_rii, pgrass$avg_PET)
# r = 0.38, p = 0.008

cor.test(pgrass$pcomm_height_rii, pgrass$avg_trans_0_200cm)
# r = -0.35, p = 0.014
cor.test(pgrass$pcomm_height_rii, pgrass$avg_trans_0_30cm)
# r = -0.36, p = 0.011
cor.test(pgrass$pcomm_height_rii, pgrass$avg_trans_30_200cm)
# r = -0.30, p = 0.037
cor.test(pgrass$pcomm_height_rii, pgrass$prop_shallow)
# r = -0.13, p = 0.392


################################################################################
# Step 6: Final figure showing RII_cover ~ amount + depth distribution of soil water

# Create binary variable for wet vs dry sites
mean(pgrass$avg_trans_0_200cm)
pgrass$dry <- as.numeric(pgrass$avg_trans_0_200cm < mean(pgrass$avg_trans_0_200cm))

# Create a second binary variable for wet vs. dry sites based on median
pgrass$dry_median <- as.numeric(pgrass$avg_trans_0_200cm < median(pgrass$avg_trans_0_200cm))

# Look at pcomm_cover_rii ~ prop_shallow + dry
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow*pgrass$dry))
# R2 = 0.32, p < 0.001, interaction not significant
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$dry))
# R2 = 0.32, p < 0.001

# Now try using median of avg_trans
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow*pgrass$dry_median))
# R2 = 0.31, p = 0.001, interaction not significant
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$dry_median))
# R2 = 0.31, p < 0.001

# Now look using average trans
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow*pgrass$avg_trans_0_200cm))
# R2 = 0.37, p < 0.001, interaction not significant
summary(lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$avg_trans_0_200cm))
# R2 = 0.37, p < 0.001

# Now with linear models only
m <- lm(pgrass$pcomm_cover_rii~pgrass$prop_shallow+pgrass$dry)

png(file.path(figdir,"Fig5b.png"), width = 4, height = 4, units = "in", res = 600)
par(mar= c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(pgrass$pcomm_cover_rii~pgrass$prop_shallow, pch = 16, 
     col = c("#01665eb3", "#8c510ab3")[pgrass$dry+1],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(-1,1),
     ylab =  "RII: Competitive (-) to Facilitative (+)", 
     xlab = "")
#points(cover1$dominance~plots$prop_shallow, cex = 2, lwd = 2)
mtext("Proportion of total transpiration from\nshallow soils layers (0-30 cm)", 1, line = 2)
abline(h=0, lty = 2)
curve(m$coefficients[1] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#01665e") # wet plots
curve(m$coefficients[1] + m$coefficients[3] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#8c510a") # dry plots
mtext("             Dry: Y = 1.91 - 3.55*X", side = 1, line = -2.2, adj = 0)
mtext("             Wet: Y = 1.67 - 3.55*X", side = 1, line = -1, adj = 0)
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

png(file.path(figdir,"Fig5a.png"), width = 4, height = 4, units = "in", res = 600)
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
#mtext("  r = -0.2, p = 0.174", side = 1, line = -1.2, adj = 0)
text("r = -0.2   ", x = 650, y = 0.97)
text("p = 0.174", x = 650, y = 0.84)
dev.off()

