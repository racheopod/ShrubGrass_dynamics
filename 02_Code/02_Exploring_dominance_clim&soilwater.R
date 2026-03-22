################################################################################
# Purpose: Evaluate plant functional type dominance across 51 shrub-grass interaction 
#          sites in sagebrush habitat across the western U.S.
#
#
# Rachel R. Renne
# March 26, 2025
# Updated: July 10, 2025
################################################################################

# Load libraries
library(dplyr)
library(terra)

# Set up directories
datadir <- '01_Data'
codedir <- '02_Code'
figdir <- '03_Figures'

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
codes1 <- codes[,c("code","GrowthForm")]
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
# Step 2: Calculate monthly PPT, TEMP correlation for Sala et al. 1997 model

# Check plots vs cover1 order
table(plots$plot == cover1$plot)

# Now calculate correlation
cover1$CORR_PT <- NA
for (i in 1:51){
  m <- cor.test(t(plots[i,24:35]), t(plots[i,36:47]))
  cover1$CORR_PT[i] <- m$estimate
}

# Calculate water holding capacity (theta_S) using the following equations from:
# Saxton, K.E., Rawls, W.J., 2006. Soil Water Characteristic Estimates by Texture 
# and Organic Matter for Hydrologic Solutions. Soil Science Society of America 
# Journal 70, 1569. https://doi.org/10.2136/sssaj2005.0117

# theta_S = theta_33 + theta_S_33 - 0.097S + 0.043
# theta_33 = theta_33t + (1.283*(theta_33t^2) 0.374*theta_33t - 0.015)
# theta_33t = -0.251S + 0.195C + 0.011OM + 
#             0.006*(S*OM) - 0.027*(C*OM) +
#             0.452*(S*C) + 0.299
# theta_S_33 = theta_S_33t + (0.636*theta_S_33t - 0.107)
# theta_S_33t = 0.278*S + 0.034*C + 0.022*OM - 
#               0.018*(S*OM) - 0.027*(C*OM) -
#               0.584*(S*C) + 0.078


# Set OM to 0.0025 (most dryland soils have low OM)
OM = 0.0025
theta_33t = -0.251*plots$sand + 0.195*plots$clay + 0.011*OM + 
             0.006*(plots$sand*OM) - 0.027*(plots$clay*OM) +
             0.452*(plots$sand*plots$clay) + 0.299

theta_33 = theta_33t + (1.283*(theta_33t^2) - 0.374*theta_33t - 0.015)

theta_S_33t = 0.278*plots$sand + 0.034*plots$clay + 0.022*OM - 
              0.018*(plots$sand*OM) - 0.027*(plots$clay*OM) -
              0.584*(plots$sand*plots$clay) + 0.078

theta_S_33 = theta_S_33t + (0.636*theta_S_33t - 0.107)

theta_S = theta_33 + theta_S_33 - 0.097*plots$sand + 0.043

# Add Water holding capacity (WHC) to cover
table(cover1$plot == plots$plot)
cover1$WHC <- theta_33

# Set up color scheme to display dominance
cols <- rev(c("#40004bb3","#40004bb3","#762a83b3","#9970abb3","#c2a5cfb3","#808080B3",
              "#a6dba0b3","#5aae61b3","#1b7837b3","#00441bb3"))

# Make cuts in dominance for display
cover1$dom_bin <- cut(cover1$dominance, breaks = c(seq(-1,0,by = 0.25),0.0000001,seq(0,1,by=0.25)[-1],1.1))

# Test correlations between dominance and climate & soil variables
cor.test(cover1$dominance, cover1$CORR_PT)
# r = -0.36, p = 0.009
cor.test(cover1$dominance, cover1$WHC)
# r = -0.08, p = 0.595
cor.test(cover1$dominance, plots$sand)
# r = 0.13, p = 0.351
cor.test(cover1$dominance, plots$clay)
# r = -0.024, p = 0.870
cor.test(cover1$dominance, plots$map)
# r = -0.16, p = 0.266
cor.test(cover1$dominance, plots$mat)
# r = 0.422, p = 0.002

# Model dominance ~ WHC + CORR_PT (based on Sala et al. 1997 model)
summary(lm(cover1$dominance~cover1$WHC+cover1$CORR_PT))
# R2 = 0.13, p = 0.033

# Model dominance ~ WHC * CORR_PT (based on Sala et al. 1997 model)
summary(lm(cover1$dominance~cover1$WHC*cover1$CORR_PT))
# R2 = 0.14, p = 0.073

# Create figure 4a and save to file
png(file.path(figdir,"Fig4a.png"), width = 4, height = 4, units = "in", res = 600)
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

# Make a map of dominance across the West
# Load sagebrush region from:
# Renne, R.R., Schlaepfer, D.R., Palmquist, K.A., Lauenroth, W.K., Bradford, J.B., 
# 2024. Estimating multivariate ecological variables at high spatial resolution 
# using a cost‐effective matching algorithm. Ecosphere 15, e4811. 
# https://doi.org/10.1002/ecs2.4811
art <- rast(file.path(datadir, 'CoreARTR_combined_DayMet_cropped&trimmed.tif'))

# Make file
png(file.path(figdir,"Fig3a.png"), width = 4, height = 4, units = "in", res = 600)
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

# Create legend for color scheme of dominance for map
png(file.path(figdir,"Fig3a_legend_v1.png"), width = 4, height = 4, units = "in", res = 600)
par(mar = c(1,1,1,1), las = 0, mgp = c(1,0.1,0))
plot(c(1:10)~1, pch = 16, col = "white", yaxt = "n", xaxt = "n", bty = "n",
     xlab = "", ylab ="")
for (i in 1:(length(cols)-1)){
  polygon(c(5,5,5.5,5.5), c(i,i+1,i+1,i), col = cols[i], lwd = 2)
}
axis(side = 2, at = c(2,9), labels = c("Grass","Shrub"), pos = 5, tick = F,
     cex.axis = 1.7)
par(mgp = c(1,0.5,0), tcl = -0.2)
axis(side = 4, at = seq(1.5,9.5), labels = c("-1","-0.75","-0.5","-0.25","0",
                                             "0.25","0.5","0.75","1"), pos = 5.53, tick = T,
     cex.axis = 1.2)
dev.off()

# Show points in climate space
png(file.path(figdir,"Fig3b.png"), width = 4, height = 4, units = "in", res = 600)
par(mar = c(2,2,1,1), tcl = 0.1, mgp = c(1,0.1,0))
plot(plots$map~plots$mat, pch = 16, cex=2,
       col = cols[cover1$dom_bin],
     xlab = "Mean annual temperature (C)",
     ylab = "Mean annual precipitation (mm)")
points(plots$map~plots$mat, pch = 1, cex=2, lwd = 2)
dev.off()

# Check plot order in cover1 and plots dfs
table(plots$plot == cover1$plot)

# Make model predicting dominance from prop_shallow and total transpiration
summary(lm(cover1$dominance ~ plots$prop_shallow*plots$avg_trans_0_200cm))
# R2 = 0.30, p = 0.001, interaction not significant

summary(lm(cover1$dominance ~ plots$prop_shallow+plots$avg_trans_0_200cm))
# R2 = 0.30, p = 0.0002, avg_trans_0_200cm not significant


################################################################################
# Step 3: Check for dung effects

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

#Reorder
dung2 <- dung2a[order(dung2a$plot),]

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
cor.test(plots$avg_trans_0_200cm, dung2$cow)
# r = 0.3130, p = 0.0253
cor.test(plots$avg_trans_0_200cm, dung2$horse)
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
# Step 4: Correlations between dominance and 1) climate, 2) transpiration, 3) vegetation

# Check that plots are aligned
table(cover1$plot == plots$plot)

cor.test(cover1$dominance, plots$map)
# r = -0.16, p = 0.266
cor.test(cover1$dominance, plots$mat)
# r = 0.42, p = 0.002
cor.test(cover1$dominance, plots$avg_PET)
# r = 0.39, p = 0.005

cor.test(cover1$dominance, plots$avg_trans_0_200cm)
# r = -0.26, p = 0.061
cor.test(cover1$dominance, plots$avg_trans_0_30cm)
# r = -0.42, p = 0.002
cor.test(cover1$dominance, plots$avg_trans_30_200cm)
# r = -0.09, p = 0.532
cor.test(cover1$dominance, plots$prop_shallow)
# r = -0.45, p = 0.001

cor.test(cover1$dominance, cover1$SH)
# r = 0.18, p = 0.198
cor.test(cover1$dominance, cover1$PG)
# r = -0.90, p < 0.001

################################################################################
# Step 5: Final figure showing dominance ~ amount + depth distribution of soil water

# Create binary variable for wet vs dry sites
mean(plots$avg_trans_0_200cm)
# 146
median(plots$avg_trans_0_200cm)
# 134
plots$dry <- as.numeric(plots$avg_trans_0_200cm < mean(plots$avg_trans_0_200cm))

# Create a second binary variable for wet vs. dry sites based on median
plots$dry_median <- as.numeric(plots$avg_trans_0_200cm < median(plots$avg_trans_0_200cm))

# Look at dominance ~ prop_shallow + dry
summary(lm(cover1$dominance~plots$prop_shallow*plots$dry))
# R2 = 0.29, p = 0.0008, interaction not significant
summary(lm(cover1$dominance~plots$prop_shallow+plots$dry))
# R2 = 0.29, p < 0.001

# Now try using median of avg_trans
summary(lm(cover1$dominance~plots$prop_shallow*plots$dry_median))
# R2 = 0.28, p = 0.001, interaction not significant
summary(lm(cover1$dominance~plots$prop_shallow+plots$dry_median))
# R2 = 0.27, p < 0.001

# Now look using average trans
summary(lm(cover1$dominance~plots$prop_shallow*plots$avg_trans_0_200cm))
# R2 = 0.30, p < 0.001, interaction not significant
summary(lm(cover1$dominance~plots$prop_shallow+plots$avg_trans_0_200cm))
# R2 = 0.30, p < 0.001

# Now with linear models only
m <- lm(cover1$dominance~plots$prop_shallow+plots$dry)

png(file.path(figdir,"Fig4b.png"), width = 4, height = 4, units = "in", res = 600)
par(mar= c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1, mfrow = c(1,1))
plot(cover1$dominance~plots$prop_shallow, pch = 16, 
     col = c("#01665eb3", "#8c510ab3")[plots$dry+1],
     cex = 2,
     xlim = c(0.3,0.7), ylim = c(-1,1),
     ylab =  "Grass (-) to shrub (+) dominance", 
     xlab = "")
#points(cover1$dominance~plots$prop_shallow, cex = 2, lwd = 2)
mtext("Proportion of total transpiration from\nshallow soil layers (0-30 cm)", 1, line = 2)
abline(h=0, lty = 2)
curve(m$coefficients[1] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#01665e") # wet plots
curve(m$coefficients[1] + m$coefficients[3] + m$coefficients[2]*x, 0.3, 0.7, add = TRUE,
      lwd = 3, col = "#8c510a") # dry plots
mtext("             Low moisture: Y = 1.65 - 2.56*X", side = 1, line = -2.2, adj = 0)
mtext("             High moisture: Y = 1.39 - 2.56*X", side = 1, line = -1, adj = 0)
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

