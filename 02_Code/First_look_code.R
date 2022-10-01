# Make some preliminary figures with data from summer:

neighbors <- read.csv("C:/Users/Rachel Renne/Dropbox/Doctoral Projects/Data/Entered Data/neighbors_FINAL.csv")
head(neighbors)
unique(neighbors$species)
table(neighbors$species)

# Calculate average distance between target and neighbor shrubs for each plot:
# Make a unique identify for each target shrub:
tshrub <- paste0(neighbors$plot,"_",neighbors$transect,"_",neighbors$point)
avgdist <- tapply(neighbors$distance, INDEX = tshrub, FUN = mean)
par(mar = c(2.5,2.5,1,1), mgp = c(1.2,0.1,0), tcl = 0.2, las = 1)
hist(avgdist, breaks = 50, col = "white",
     xlab = "Mean distance between target and neighbor shrubs (cm)",
     main = "Neighbor distances")
abline(h = 0)
box()

# Now average over the plot:
tshrub1 <- gsub("_.*","",names(avgdist))
plotavgdist <- tapply(avgdist, INDEX = tshrub1, FUN = mean)

# Read in coordinates
coords <- read.csv("C:/Users/Rachel Renne/Dropbox/Doctoral Projects/Data/Entered Data/coordinates.csv")
coords$plot <- gsub("-P.*","", coords$ident)
coords$plot <- gsub("-C.*","", coords$plot)
coords$plot <- gsub("PC.*","", coords$plot)

# Make dataframe with coords, avgdist, plot:
dat <- data.frame(plot = toupper(names(plotavgdist)), plotavgdist = plotavgdist)
dat1 <- merge(dat, coords[,c(7,2,3,6)], by = "plot")

par(mar = c(0,2,0,2))
plot(dat1$Latitude ~ dat1$Longitude, cex = dat1$plotavgdist*0.01,
     xlim = c(-125, -103), ylim = c(32,50), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "", ylab = "", pch = 16)
maps::map("state", interior = T, add = T)

# Get climate data from daymet
library(daymetr)

# Make day/month list
month <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))

# Make dataframes:
ppt = data.frame(plot = dat1$plot, Jan = NA, Feb = NA, Mar = NA, Apr = NA, May = NA, Jun = NA,
                 Jul = NA, Aug = NA, Sep = NA, Oct = NA, Nov = NA, Dec = NA)
tmax = data.frame(plot = dat1$plot, Jan = NA, Feb = NA, Mar = NA, Apr = NA, May = NA, Jun = NA,
                 Jul = NA, Aug = NA, Sep = NA, Oct = NA, Nov = NA, Dec = NA)
tmin = data.frame(plot = dat1$plot, Jan = NA, Feb = NA, Mar = NA, Apr = NA, May = NA, Jun = NA,
                 Jul = NA, Aug = NA, Sep = NA, Oct = NA, Nov = NA, Dec = NA)
tabsmin = data.frame(plot = dat1$plot, Jan = NA, Feb = NA, Mar = NA, Apr = NA, May = NA, Jun = NA,
                      Jul = NA, Aug = NA, Sep = NA, Oct = NA, Nov = NA, Dec = NA)

# Run through points and make png of WL diagram:
for (i in 1:nrow(dat1)){
  dmt <- download_daymet(site = "Daymet", lat = dat1$Latitude[i], lon = dat1$Longitude[i],
                         start = 1981,
                         end = 2010,
                         internal = TRUE,
                         silent = FALSE) 
  
  # Add month column
  dmt$data$month <- month 
  
  # Create dat
  ppt[i,2:13] = tapply(dmt$data$prcp..mm.day., dmt$data$month, FUN = function(x){sum(x)/30})
  tmax[i,2:13] = tapply(dmt$data$tmax..deg.c., dmt$data$month, mean)
  tmin[i,2:13] = tapply(dmt$data$tmin..deg.c., dmt$data$month, mean)
  tabsmin[i,2:13] = tapply(dmt$data$tmin..deg.c., dmt$data$month, min)

}


# Put map and mat onto dat1
dat1$map <- apply(ppt[2:13], 1, mean)*10
dat1$mat <- apply(cbind(tmin[,2:13],tmax[,2:13]), 1, mean)

par(mar = c(3,3,1,1), tcl = 0.2, mgp = c(1.5,0.2,0), las = 1)
plot(x = dat1$map, y = dat1$plotavgdist,
     xlab = "Mean annual precipitation",
     ylab = "Average distance between shrubs")
plot(x = dat1$map, y = dat1$mat)
plot(x = dat1$mat, y = dat1$plotavgdist)
cor(dat1$plotavgdist, dat1$map)

# Read in hummocks
hummocks <- read.csv("C:/Users/Rachel Renne/Dropbox/Doctoral Projects/Data/Entered Data/hummocks_FINAL.csv")

# deal with <1 litter measurement
hummocks[hummocks$litter == "<1",4] <- 0
hummocks$litter <- as.numeric(hummocks$litter)

# Summarize mineral soil and litter by plot
litter <- tapply(hummocks$litter,hummocks$plot, mean)
mineralsoil <- tapply(hummocks$mineralsoil, hummocks$plot, mean)

# put these onto dat1
dat1$litter <- litter
dat1$soil <- mineralsoil

plot(dat1$litter~dat1$soil)
plot(dat1$litter~dat1$map)
plot(dat1$soil~dat1$map)
cor.test(dat1$soil, dat1$map)

# Bring in target shrub info
target <- read.csv("C:/Users/Rachel Renne/Dropbox/Doctoral Projects/Data/Entered Data/targetshrubs_FINAL.csv")

# calculate average diameter
target$avgdiam <- apply(target[,5:7], 1, mean)
avgdiam <- tapply(target$avgdiam, target$plot, mean)
avgpctlive <- tapply(target$percentlive, target$plot, mean)
subspecies <- tapply(target$targetshrub, target$plot, function(x){names(which.max(table(x)))})
                    
# add onto dat1
dat1$avgdiam <- avgdiam
dat1$avgpctlive <- avgpctlive
dat1$subspecies <- subspecies

plot(dat1$avgdiam~dat1$map)
plot(dat1$avgpctlive~dat1$map)
plot(dat1$plotavgdist~dat1$avgdiam, col = "white")
text(y = dat1$plotavgdist, x = dat1$avgdiam, labels = subspecies, cex = 0.7, offset = 0)

plot(dat1$litter~dat1$soil, col = "white")
text(y = dat1$litter, x = dat1$soil, labels = subspecies, cex = 0.7, offset = 0)
plot(dat1$litter~dat1$map)
plot(dat1$soil~dat1$map, col = "white")
text(y = dat1$soil, x = dat1$map, labels = subspecies, cex = 0.7, offset = 0)

plot(dat1$soil~dat1$mat, col = "white")
text(y = dat1$soil, x = dat1$mat, labels = subspecies, cex = 0.7, offset = 0)

plot(dat1$soil~dat1$altitude)


# Some "fun" plots
plot(dat1$plotavgdist~dat1$avgdiam,
     xlab = "Mean target shrub diameter (cm)",
     ylab = "Mean distance between shrubs (cm)")
plot(dat1$plotavgdist~dat1$map,
     xlab = "Mean annual precipitation (mm)",
     ylab = "Mean distance between shrubs (cm)")


plot(dat1$soil~dat1$map,
     xlab = "Mean annual precipitation (mm)",
     ylab = "Mean depth of hummock soil (cm)")
plot(dat1$avgdiam~dat1$map,
     ylab = "Mean target shrub diameter (cm)",
     xlab = "Mean annual precipitation")

plot(dat1$soil~dat1$plotavgdist,
     xlab = "Mean distance between shrubs (cm)",
     ylab = "Mean depth of hummock soil (cm)")
plot(dat1$soil~dat1$avgdiam,
     xlab = "Mean target shrub diameter (cm)",
     ylab = "Mean depth of hummock soil (cm)")
