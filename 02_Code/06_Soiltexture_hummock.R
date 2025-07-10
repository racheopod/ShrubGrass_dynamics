dat <- read.csv("C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2021/01_Raw_data/soiltexture_2021_preliminarycalcs.csv")
head(dat)

dat$layer <- paste0(dat$plot, "_layer", round(dat$depth_0/10)+1)

for (i in 1:nrow(dat)){
  if (grepl("NA", dat$layer[i])){
    dat$layer[i] <- paste0(dat$plot[i],"_","hummock")
  }
}

dat2 <- dat[,c("plot","sand","clay","silt","layer")]

dat2 <- dat2[grepl("layer1|hummock", dat2$layer), ]
dat2 <- dat2[order(dat2$plot),]
table(dat2$plot)

# Remove RkyMtn5-015 bc it didn't have hummocks
dat2 <- dat2[dat2$plot != "RkyMtn5-015",]

hum <- dat2[grepl("hummock",dat2$layer),]
layer1 <- dat2[grepl("layer1",dat2$layer),]

for (i in 1:nrow(hum)){
  hum$sand_diff[i] <- hum$sand[i] - layer1[layer1$plot == hum$plot[i], ]$sand
  hum$clay_diff[i] <- hum$clay[i] - layer1[layer1$plot == hum$plot[i], ]$clay
  hum$silt_diff[i] <- hum$silt[i] - layer1[layer1$plot == hum$plot[i], ]$silt
}


library(MASS)
par(mfrow=c(3,1), mgp=c(1,0.1,0),tcl=0.1, mar = c(2,2,1,1))
truehist(hum$sand_diff, xlim = c(-0.45, 0.45),
         xlab = "Hummock - Interspace sand content",
         breaks = seq(-0.45,0.45,by = 0.05), col= "grey")
abline(h=0)
box()
truehist(hum$clay_diff, xlim = c(-0.45, 0.45),
         xlab = "Hummock - Interspace clay content",
         breaks = seq(-0.45,0.45,by = 0.05), col= "grey")
abline(h=0)
box()
truehist(hum$silt_diff, xlim = c(-0.45, 0.45),
         xlab = "Hummock - Interspace silt content",
         breaks = seq(-0.45,0.45,by = 0.05), col= "grey")
abline(h=0)
box()

##########################################################
# Absolute content:

# Add sand grams
dat$Silt_g <- 50 - dat$Sand_g - dat$Clay_g

dat3 <- dat[,c("plot","Sand_g","Clay_g","Silt_g","layer")]

dat2 <- dat2[grepl("layer1|hummock", dat2$layer), ]
dat2 <- dat2[order(dat2$plot),]
table(dat2$plot)

# Remove RkyMtn5-015 bc it didn't have hummocks
dat2 <- dat2[dat2$plot != "RkyMtn5-015",]

hum <- dat2[grepl("hummock",dat2$layer),]
layer1 <- dat2[grepl("layer1",dat2$layer),]

for (i in 1:nrow(hum)){
  hum$sand_diff[i] <- hum$sand[i] - layer1[layer1$plot == hum$plot[i], ]$sand
  hum$clay_diff[i] <- hum$clay[i] - layer1[layer1$plot == hum$plot[i], ]$clay
  hum$silt_diff[i] <- hum$silt[i] - layer1[layer1$plot == hum$plot[i], ]$silt
}


library(MASS)
par(mfrow=c(3,1), mgp=c(1,0.1,0),tcl=0.1, mar = c(2,2,1,1))
truehist(hum$sand_diff, xlim = c(-0.45, 0.45),
         xlab = "Hummock - Interspace sand content",
         breaks = seq(-0.45,0.45,by = 0.05), col= "grey")
abline(h=0)
box()
truehist(hum$clay_diff, xlim = c(-0.45, 0.45),
         xlab = "Hummock - Interspace clay content",
         breaks = seq(-0.45,0.45,by = 0.05), col= "grey")
abline(h=0)
box()
truehist(hum$silt_diff, xlim = c(-0.45, 0.45),
         xlab = "Hummock - Interspace silt content",
         breaks = seq(-0.45,0.45,by = 0.05), col= "grey")
abline(h=0)
box()
