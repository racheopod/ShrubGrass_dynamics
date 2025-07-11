################################################################################
# Purpose: Vegetation sensitivity analysis for RII ~ prop_shallow results
# 
# Rachel R. Renne
# January 17, 2025
# Updated: July 11, 2025
################################################################################

# Approach: Results for each of the 51 plots with cover settings from the other 50 plots.
# Part 1: Randomly select results for each plot 10,000 times and calculate corr.coef
# between prop_shallow and RII. Compare to site-specific vegetation results
# Part 2: Use same vegetation settings for all 51 plots (51 different settings,
# one for each plot). Calculate correlation between prop_shallow and RII.
# Compare to site-specific vegetation results.

# Load relevant libraries
library("tidyr")
library(Ternary)

# Set up directories
datdir <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/01_Data/Data for submission"
figdir <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/03_Figures/Figures for submission"
codedir <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/02_Code"

################################################################################
# STEP 1: Get simulation results with site-specific vegetation cover settings

plots <- read.csv(file.path(datdir, "plots.csv"))

################################################################################
# STEP 2: Calculate RIIs from plant community data

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

# Merge with plots
results1 <- merge(plots, pg[,c(1,17)], by = "plot", all.x = TRUE)
 
# Look at figure of RII ~ prop_shallow
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(results1$pcomm_cover_rii~results1$prop_shallow, 
     xlim = c(0.3, 0.7), ylim = c(-1,1),
     col = rgb(0.5,0.5,0.5,0.7), pch = 16, cex = 2,
     xlab = "Proportion shallow transpiration (0-30 cm)",
     ylab = "RII (competition to facilitation)")
abline(h=0, lty = 2)

# Save results (correlation)
exp_results <- cor.test(results1$pcomm_cover_rii, results1$prop_shallow)

################################################################################
# STEP 3: Get vegetation cover sensitivity simulation results

# Read in simulation results
results_wide <- read.csv(file.path(datdir,"Vegetation_sensitivity_simulation_results.csv"))

# Merge with results1 (don't keep veg info here)
results2 <- merge(results1, results_wide[,c(1,8:58)], by = "plot")

# How often does range of prop_shallow not include experimental prop_shallow?
res1 <- data.frame(plot = results2$plot, prop_shallow = results2$prop_shallow,
                   min = apply(results2[,49:99],1,min), max = apply(results2[,49:99],1,max)) 

res1$inrange <- as.numeric(round(res1$prop_shallow - res1$min,4) >= 0 & round(res1$max - res1$prop_shallow,4) >= 0)
sum(res1$inrange)
# 51 --> All capture the actual simulated value--though it is sometimes the most extreme value

# Create a boxplots of prop_shallow for each plot
# Custom boxplot function
placebox <- function(y, x, width = 0.02, lwd = 1.5, col = "black",
                     border = border){
  # Calculate necessary stats
  sts = boxplot.stats(x)
  
  # Make box
  polygon(x = c(sts$stats[2], sts$stats[2], sts$stats[4],sts$stats[4]),
          y = c((y - width), (y + width), (y + width), (y - width)),
          lwd = lwd, border = border)
  # Add median
  lines(x = c(sts$stats[3],sts$stats[3]), y = c((y - width),(y + width)), 
        lwd = lwd, col = col)
  
  # add wiskers
  arrows(x0 = sts$stats[2], x1 = sts$stats[1], y0 = y, y1 = y, angle = 90,
         lwd = lwd, length = width*2, col = col)
  arrows(x0 = sts$stats[4], x1 = sts$stats[5], y0 = y, y1 = y, angle = 90,
         lwd = lwd, length = width*2, col = col)
  # add outliers
  if (length(sts$out > 0)){
    points(x = sts$out, y = rep(y, length(sts$out)), pch = 16, 
           col = rgb(0.5,0.5,0.5,0.5), cex = 1)
  }
}

# Create figure
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(1, xlim = c(0.3,0.7), ylim = c(-1,1), xlab = "Proportion shallow transpiration",
     ylab = "RII (competitive to facilitative)", col = "white")
abline(h=0, lty = 3)
# loop through results and "placebox"
for (i in 1:51){
  y = results2$pcomm_cover_rii[i]
  x = t(results2[i,49:99])
  placebox(y, x, width = 0.01, lwd = 1.5, col = rgb(0.5,0.5,0.5,0.7),
                       border = rgb(0.5,0.5,0.5,0.7))
}

# Another function just to add whiskers
placewhiskers <- function(y, x, width = 0.02, lwd = 1.5, col = "black"){
  # Calculate necessary stats
  sts = boxplot.stats(x)
  
  # add wiskers
  arrows(x0 = sts$stats[3], x1 = sts$stats[1], y0 = y, y1 = y, angle = 90,
         lwd = lwd, length = width*2, col = col)
  arrows(x0 = sts$stats[3], x1 = sts$stats[5], y0 = y, y1 = y, angle = 90,
         lwd = lwd, length = width*2, col = col)
}

# Create another figure that adds points to the experimental results
par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(1, xlim = c(0.3,0.7), ylim = c(-1,1), xlab = "Proportion shallow transpiration",
     ylab = "RII (competitive to facilitative)", col = "white")
abline(h=0, lty = 3)
# loop through results and "placewhiskers"
for (i in 1:51){
  y = results2$pcomm_cover_rii[i]
  x = t(results2[i,49:99])
  placewhiskers(y, x, width = 0.01, lwd = 1.5, col = rgb(0.5,0.5,0.5,0.7))
}
points(results2$pcomm_cover_rii~results2$prop_shallow, pch = 16,
       cex = 1.5, col = rgb(0.5,0.5,0.5,0.5))


# Create df to store results
results <- data.frame(coeff = rep(NA,10000), pvalue = rep(NA,10000))

# Save plot that will contain all lines
png(file.path(figdir,"Chapter1_Supp4FigS1a.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(1, xlim = c(0.3,0.7), ylim = c(-1,1), xlab = "",
     ylab = "RII (competitive to facilitative)", col = "white")
mtext("Proportion of total transpiration from\nshallow soil layers (0-30 cm)", 1, line = 2)
abline(h=0, lty = 3)

# Loop through and repeat 10000 times
for (i in 1:10000){
  # Report progress
  if (i %in% seq(100,10000, by = 100)){
    print(paste0("Now working on iteration ",i,"."))
  } 
  
  # Get random selection of results and calculate correlation
  x1 <- sample(1:51, 51, TRUE)
  x2 <- apply(cbind(results2[,49:99],x1),1,function(x){x[x[52]]})
  
  # Calculate correlation coefficient
  cc <- cor.test(results2$pcomm_cover_rii, x2)
  
  # Add results to dataframe
  results$coeff[i] <- cc$estimate
  results$pvalue[i] <- cc$p.value
  
  # Now make LM
  m <- lm(results2$pcomm_cover_rii~x2)
  curve(m$coefficients[1]+x*m$coefficients[2], from = 0.3, to = 0.7, 
         col = c(rgb(0.5,0.5,0.5,0.01),rgb(1,0.2,0,0.01))[(cc$p.value <= 0.05)+1], lwd = 2, add = TRUE)
}


# Now see how often results are significant
# Results vary very slightly due to random draws
prop_sig <- round(mean(results$pvalue <= 0.05)*100)
# 0.2426
# 0.2698
# 0.2539

# See how often coefficient is >= 0
# Results vary very slightly due to random draws
prop_gt0 <- min(99.9, round(100 - 100*table(results$coef >= 0)[2]/10000, 1))
# 0.0009
# 0.0007
# 0.0004

# Add to plot
legend("topright", bty = "n", legend = c(paste0("p < 0.05 for ",prop_sig,"% of samples"),
                                         paste0("r < 0 for ",prop_gt0,"% of samples")))
legend("bottomleft", bty = "n", legend = "n = 10000")
# Manually export plot
dev.off()

results[results$coeff == max(results[results$pvalue <= 0.05,]$coeff),]
# smallest magnitude correlation that is significant: 
# -0.2759
# -0.2816
# Results vary very slightly due to random draws

# What proportion of simulations have coeff < 0
mean(results$coeff < 0)


# Determine percentage results <= site-specific vegetation results
ests <- mean(results$coeff <= exp_results$estimate)
# 1e-04
# Results vary very slightly due to random draws


# Create histogram of correlation coefficients
png(file.path(figdir,"Chapter1_Supp4FigS1b.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
hist(results$coeff, breaks = 50,
     main = "",
     xlab = "",
     xlim = c(-0.53,0.04))
abline(v = exp_results$estimate, lty = 2, col = "red", lwd = 2)
mtext("Correlation coefficient", 1, line = 1.5)
box()
legend("topright", bty = "n", legend = c("Estimate < 100%","of samples"), 
       lty = 2, col = c("red","white"), lwd = 2, cex = 0.7)
abline(h=0)
dev.off()


################################################################################
# STEP 4: Now look at results if we set vegetation to be the same for all plots.

# Create df to store results
sameveg <- data.frame(vegplot = results2$plot, coeff = rep(NA, 51), pvalue = NA)

# Plot linear relationships when veg is the same for all plots
# Create figure 
png(file.path(figdir,"Chapter1_Supp4FigS2a.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(3,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(1, xlim = c(0.3,0.7), ylim = c(-1,1), xlab = "",
     ylab = "RII (competitive to facilitative)", col = "white")
mtext("Proportion of total transpiration from\nshallow soil layers (0-30 cm)", 1, line = 2)
abline(h=0, lty = 3)


# Loop through results2 df and calculate values
for (i in 1:51){
  
  # Get results from one veg type across all plots.
  x2 <- results2[,48+i]
  
  # Calculate correlation coefficient
  cc <- cor.test(results2$pcomm_cover_rii, x2)
  
  # Add results to dataframe
  sameveg$coeff[i] <- cc$estimate
  sameveg$pvalue[i] <- cc$p.value
  
  # Now make LM & add line to the plot
  m <- lm(results2$pcomm_cover_rii~x2)
  curve(m$coefficients[1]+x*m$coefficients[2], from = 0.3, to = 0.7, 
        col = c(rgb(0.5,0.5,0.5,0.3),rgb(1,0.2,0,0.3))[(cc$p.value <= 0.05)+1], lwd = 2, add = TRUE)
}


# Results vary due to random draws
prop_sig <- round(mean(sameveg$pvalue <= 0.05)*100)
# 0.39
# 0.412
summary(sameveg$coeff)
# -0.31 to -0.19 (Nothing as extreme as with site-specific veg parameters)
# -0.32 to -0.20
prop_lt0 <- round(mean(sameveg$coeff < 0)*100)


legend('bottomleft',bty="n",legend=c(expression(p > 0.05), expression(p <= 0.05)),
       col = c(rgb(0.5,0.5,0.5,0.7),rgb(1,0.2,0,0.7)), lwd = 2)
legend('topright',bty="n",legend=c(paste0("p < 0.05 for ",prop_sig,"% of samples"),
                                   paste0("r < 0 for ",prop_lt0,"% of samples")))
dev.off()


hist(sameveg$coeff, breaks = 20)
abline(v = max(sameveg[sameveg$pvalue <= 0.05,]$coeff), lty = 2, lwd = 2, col = "red")
abline(v = exp_results$estimate, lty = 2, lwd = 2, col = "blue")

sameveg[which.min(sameveg$coeff),]
sameveg[which.max(sameveg$coeff),]

# Look at which ones were significant
par(las = 1, mar = c(2,6,1,1))
boxplot(results2[,49:99], horizontal = TRUE, cex.axis = 0.5, 
        col = c("grey", "red")[(sameveg$pvalue <= 0.05) + 1])


################################################################################
# STEP 5: Explore which vegetation settings resulted in significant correlations

# Make matrix of results_wide
veg1 <- data.frame(shrub = results_wide$shrubs, 
                   herb = c(results_wide$c3pgrass+results_wide$c4pgrass+results_wide$agrass+results_wide$pforbs),
                   bareground = results_wide$bareground)

# Save ternary plot
png(file.path(figdir,"Chapter1_Supp4FigS2b.png"), width = 4, height = 4, units = "in", res = 300)
par(mar = c(0,0,0,0))
TernaryPlot(alab = "Shrub", blab = "Herbaceous", clab = "Bareground")
TernaryPoints(as.matrix(veg1), pch = 16, 
              col = c("grey", "red")[(sameveg$pvalue <= 0.05) + 1])
legend(x = 0.22, y = 0.88,bty="n",legend=c(expression(p > 0.05), expression(p <= 0.05)),
       col = c(rgb(0.5,0.5,0.5,0.7),rgb(1,0.2,0,0.7)), pch = 16)
dev.off()
