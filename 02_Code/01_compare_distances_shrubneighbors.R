
##################################### 
# Purpose: Check distances between target and neighbor shrubs measured by two different people
#
#
#
# Rachel R. Renne
# September 14, 2022
#####################################

datadir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data'

# Read in distances
dists <- read.csv(file.path(datadir,"distances.csv"))
# rename "distance" column
names(dists)[4] <- "dist_rrr"

# Read in neighbors
neighbors <- read.csv(file.path(datadir, "neighbors.csv"))

# Pull out the neighbors that were used for plant community measurements
# Keep only relevant columns
neighbors1 <- neighbors[neighbors$neighbor == 1, c(1,2,3,7)]

# Merge data frames
dists_merged <- merge(dists, neighbors1, by = c("plot","transect","point"))

# Read in plots
plots <- read.csv(file.path(datadir, "plots.csv"))
# limit to plotname and date
plots1 <- plots[,c(1,5)]

# Merge plots and dists_merged
dists1 <- merge(plots1, dists_merged, by = "plot")

# Make dataframe showing each technician x date
techs <- data.frame(date = sort(unique(dists1$date)), tech = NA)
techs$tech[1:7] <- "AR"
techs$tech[8:14] <- "PF"
techs$tech[15:38] <- "JG"

# Merge techs with dists1
dists2 <- merge(dists1, techs, by = "date")

##############################################
# 1. Compare across all plots, all techs

# First, plot
plot(dists2$distance~dists2$dist_rrr)

# Model with all techs
m <- lm(dists2$distance ~ dists2$dist_rrr)
summary(m)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      8.29063    1.28520   6.451 2.24e-10 ***
#   dists2$dist_rrr  0.90216    0.01008  89.476  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 13.01 on 621 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.928,	Adjusted R-squared:  0.9279 
# F-statistic:  8006 on 1 and 621 DF,  p-value: < 2.2e-16

# Use model to predict to show relationship
x <- seq(42,389,1)
y <- m$coefficients[1] + m$coefficients[2]*x

# Plot
par(mar = c(2,2,2,1), mgp = c(1,0.2,0), tcl = 0.1)
plot(dists2$distance~dists2$dist_rrr,
     main = "Distance measured by two observers",
     xlab = "Distance from plant community protocol",
     ylab = "Distance from shrub neighbor protocol",
     xlim = c(40,400), ylim = c(40,400),
     pch = 16, lwd = 2, col = rgb(0.5,0.5,0.5,0.5))#col = c("red","blue","black")[as.factor(dists2$tech)])
lines(y~x, lwd = 2)
abline(0,1, lty = 2)

# Model with tech as variable
m1 <- lm(dists2$distance ~ dists2$dist_rrr * dists2$tech)
summary(m1)
 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    28.20225    3.50459   8.047 4.36e-15 ***
#   dists2$dist_rrr                 0.68849    0.03304  20.837  < 2e-16 ***
#   dists2$techJG                 -21.59980    3.82162  -5.652 2.42e-08 ***
#   dists2$techPF                 -28.25210    4.87854  -5.791 1.11e-08 ***
#   dists2$dist_rrr:dists2$techJG   0.22493    0.03483   6.458 2.15e-10 ***
#   dists2$dist_rrr:dists2$techPF   0.31285    0.04478   6.986 7.33e-12 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 12.46 on 617 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9344,	Adjusted R-squared:  0.9338 
# F-statistic:  1757 on 5 and 617 DF,  p-value: < 2.2e-16

xAR <- (min(dists2[dists2$tech == "AR",]$dist_rrr,na.omit=T):max(dists2[dists2$tech == "AR",]$dist_rrr,na.omit=T))
xJG <- (min(na.omit(dists2[dists2$tech == "JG",]$dist_rrr)):max(na.omit(dists2[dists2$tech == "JG",]$dist_rrr)))
xPF <- (min(na.omit(dists2[dists2$tech == "PF",]$dist_rrr)):max(na.omit(dists2[dists2$tech == "PF",]$dist_rrr)))
AR <- m1$coefficients[1] + m1$coefficients[2]*xAR
JG <- m1$coefficients[1] + m1$coefficients[3] + (m1$coefficients[2]+m1$coefficients[5])*xJG
PF <- m1$coefficients[1] + m1$coefficients[4] + (m1$coefficients[2]+m1$coefficients[6])*xPF

# Plot
par(mar = c(2,2,2,1), mgp = c(1,0.2,0), tcl = 0.1)
plot(distance~dist_rrr, data = dists2[dists2$tech == "AR",],
     main = "Distance measured by two observers",
     xlab = "Distance from plant community protocol",
     ylab = "Distance from shrub neighbor protocol (AR)",
     xlim = c(40,400), ylim = c(40,400),
     pch = 16, lwd = 2, col = rgb(0.5,0.5,0.5,0.5))#col = c("red","blue","black")[as.factor(dists2$tech)])
lines(AR~xAR, lwd = 2)
abline(0,1, lty = 2)

par(mar = c(2,2,2,1), mgp = c(1,0.2,0), tcl = 0.1)
plot(distance~dist_rrr, data = dists2[dists2$tech == "JG",],
     main = "Distance measured by two observers",
     xlab = "Distance from plant community protocol",
     ylab = "Distance from shrub neighbor protocol (JG)",
     xlim = c(40,400), ylim = c(40,400),
     pch = 16, lwd = 2, col = rgb(0.5,0.5,0.5,0.5))#col = c("red","blue","black")[as.factor(dists2$tech)])
lines(JG~xJG, lwd = 2)
abline(0,1, lty = 2)

par(mar = c(2,2,2,1), mgp = c(1,0.2,0), tcl = 0.1)
plot(distance~dist_rrr, data = dists2[dists2$tech == "PF",],
     main = "Distance measured by two observers",
     xlab = "Distance from plant community protocol",
     ylab = "Distance from shrub neighbor protocol (PF)",
     xlim = c(40,400), ylim = c(40,400),
     pch = 16, lwd = 2, col = rgb(0.5,0.5,0.5,0.5))#col = c("red","blue","black")[as.factor(dists2$tech)])
lines(PF~xPF, lwd = 2)
abline(0,1, lty = 2)
