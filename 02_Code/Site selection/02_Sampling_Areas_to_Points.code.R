#################################
#
# Purpose: Convert Sagebrush density, climate zone, precipitation zone, and MAP to points for site selection
#
# Inputs: "CandidateSamplingAreas_2km_lowslopes_MAP-0000000000-0000000000.tif"
# 		  "CandidateSamplingAreas_2km_lowslopes_MAP-0000000000-0000023296.tif"
#		  "CandidateSamplingAreas_2km_lowslopes_czone.tif"
#		  "CandidateSamplingAreas_2km_lowslopes_pzone.tif"
# 		  "CandidateSamplingAreas_2km_lowslopes_sagebrush.tif"
#		  "landmgmt.gdb"
#
# Outputs: "CandidateSamplingAreas_2km_lowslopes_MAPmoasic.tif"
#		   "CandidateSampleingPoints.csv"
#		   "Points_czone_x_pzone.png"
#		   "CandidateSampleingPoints_publiclands_only.csv"
#		   "Renne2021_Points.csv"
#
# Rachel R. Renne
# April 12, 2021
############################

library(raster)

# Read in MAP rasters:
m1 <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_MAP-0000000000-0000000000.tif")
m2 <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_MAP-0000000000-0000023296.tif")

# Mosaic together
map <- mosaic(m1, m2, fun = max)
rm(m1,m2)
#writeRaster(map, "C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_MAPmoasic.tif")
#map <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_MAPmoasic.tif")
 
# Read in climate zones:
czone <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_czone.tif")
  
# Read in precipitation zone:
pzone <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_pzone.tif")
# Set NA value to 0
NAvalue(pzone) <- 0

# Read in sagebrush density
sage <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_sagebrush.tif")

# Raster to points (using map), make Spatial so we can use this to extract from the rest
pts1 <- rasterToPoints(map, spatial=TRUE)
pts1a <- data.frame(cbind(pts1@coords, pts1@data))
colnames(pts1a)[3] <- "MAP"

# Need to extract czone to points to get accurate representation:
pts2 <- extract(czone, pts1, method ="simple")

# Add onto pts1a
pts1a$czone <- pts2

# Now extract the rest of the rasters:
pts3 <- rasterToPoints(pzone,spatial=F)
colnames(pts3)[3] <- "pzone"

# Remove rasters and extraneous pt files
rm(map,pzone,czone,sage,pts1,pts2)

# Tack on pzone:
pts1a$pzone <- pts3[,3]
pts1a$sagebrush <- pts4[,3]

# Now extract sagebrush density from raster:
pts4 <- extract(sage, pts1a[,1:2], method = "simple")

# tack onto pts1a
pts1a$sagebrush <- pts4

#Remove lines where sagebrush = NA
pts1b <- pts1a[!is.na(pts4),]
# 23,724,893--23 million cells!!

# Save file
#saveRDS(pts1b, "C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSampleingPoints.csv")

# Remove extraneous stuff
rm(pts1a, pts3, sage,pts4)

# Read in land management:
land <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/landmgmt")

landx <- extract(land, pts1b[,1:2], method = "simple")

# Add onto pts1b
pts1b$landmgmt <- landx
rm(land)

# Get only those points in BLM, USFS, FWS, STATE, USBR
pts1c <- pts1b[landx == 3 |landx == 7 |landx == 8 | landx == 14 | landx == 19,]

# Get rid of NA's
pts1d <- pts1c[!is.na(pts1c$landmgmt),]

# merge land ownershp code with agency acronym for clarity
landy <- data.frame(landmgmt = c(3,7,8,14,19), Agency = c("BLM", "USFS", "FWS", "ST", "USBR"))

pts1e <- merge(pts1d, landy, by = "landmgmt", all.x = T)

# Summarize pzone by czone & make figures:

counts <- table(pts1e$czone, pts1e$pzone)

czonenames <- c("Pacific NW Mountains", "Dry Great Basin", "Great Plains","Southern Edge", "Colorado Plateau", "WY Basin & RM Valleys", 
                "Great Basin", "Rocky Mountains")

png("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/Points_czone_x_pzone.png", width = 12, height = 6, units = "in", res = 300)
par(mfrow=c(2,4))
for (i in c(3,6,8,1,2,7,5,4)){
  par(mar = c(3,3,2,1), mgp = c(2,0.5,0), tcl = -0.3)
  barplot(counts[i,], ylim = c(0,max(counts[i,])*1.1), main = czonenames[i],
          xlab = "Precipitation zone", ylab = "Frequency")
  legend("topright", legend = paste0("n = ", sum(counts[i,])), bty = "n",xjust = 0, cex = 1.2)
  box()
}
dev.off()

# Remove land management code
pts1f <- pts1e[,-1]

# Save to file:
#saveRDS(pts1f, "C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSampleingPoints_publiclands_only.csv")


# Remove extraneous data
rm(landy, pts1b, pts1c, pts1d, pts1e, counts, landx)

##################################
# Now select samples of sites:

# Just get pts with sagebrush density >109 (10% coverage per ~1 km2)
sagex <- pts1f[pts1f$sagebrush >= 109,]

# Create a stratum field that combines czone with pzone:
sagex$stratum <- paste0("c",sagex$czone,"p", sagex$pzone)

# Read in sage for crs
sage <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_sagebrush.tif")

# Make sagex into a spatial points df
sage1 <- SpatialPointsDataFrame(sagex[,1:2], data = sagex[,3:8], proj4string = crs(sage)) 

library(spsurvey)

# Set up the design
# Strata: c1p3,c1p4,c1p5,c1p6,c1p7,
#         c2p2,c2p3,c2p4,c2p5,
#         c3p2,c3p3,c3p4,c3p5,
#         c6p3,c6p4,c6p5,c6p6,c6p7,
#         c7p2,c7p3,c7p4,c7p5,c7p6,c7c7,
#         c8p2,c8p3,c8p4,c8p5,c8p6,c8p7,
#         c9p3,c9p4,c9p5,c9p6,c9p7,
#         cNAp2,cNAp7 # Include the extreme precip values for cells outside climate zones.
# 

rm(sagex, sage, pts1f)
gc()

thisdesign <- list(c1p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c1p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c1p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c1p6 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c1p7 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c2p2 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c2p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c2p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c2p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c3p2 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c3p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c3p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c3p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c6p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c6p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c6p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c6p6 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c6p7 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c7p2 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10))
                   
thisdesign1 <- list(c7p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c7p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c7p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c7p6 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c7p7 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c8p2 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c8p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c8p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c8p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c8p6 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c8p7 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c9p3 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c9p4 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c9p5 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c9p6 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   c9p7 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   cNAp2 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10),
                   cNAp7 = list(panel = c(PanelOne = 10), seltype = "Equal", over = 10))

samplinglocations = grts(design = thisdesign, DesignID = "Renne2021", SiteBegin = 1, type.frame = "finite", src.frame = "sp.object",
                         sp.object = sage1, stratum = "stratum")

samplinglocations1 = grts(design = thisdesign1, DesignID = "Renne2021", SiteBegin = 1, type.frame = "finite", src.frame = "sp.object",
                         sp.object = sage1, stratum = "stratum")

library(maps)
plot(samplinglocations$ycoord~samplinglocations$xcoord, col = samplinglocations@data$czone, pch = 1, ylim = c(32,49), xlim = c(-125, -100))
points(samplinglocations1$ycoord~samplinglocations1$xcoord, col = samplinglocations1@data$czone, pch = 2)
map('state', interior = T, add = T)


# Gather relevant data and save to file
samplepts <- data.frame(rbind(samplinglocations@data, samplinglocations1@data))

# Get pts numbered sequentially
fixSiteID <- paste0("Renne2021", "-", c(1:nrow(samplepts)))
samplepts$siteID[100:nrow(samplepts)] <- fixSiteID[100:nrow(samplepts)]

# Merge names of zones with czone:
zonekey <- data.frame(czonename = czonenames <- c("PNWMtn", "DryGB", "GP",NA, "SE", "CP", "WYBas", 
                                                  "GB", "RkyMtn"), czone = c(1:9))

samplepts1 <- merge(samplepts, zonekey, by = "czone", all.x = T)

# Now create sequential plot ID's
plotIDs <- rep(c(1:20), nrow(samplepts1)/20)

zeroseq <- rep(c(rep("00",9), rep("0",11)),nrow(samplepts1)/20)

samplepts1$plotseq = plotIDs

samplepts1$PlotID <- paste0(samplepts1$czonename, samplepts1$pzone, "-", zeroseq, samplepts1$plotseq)

# Save to file:
#write.csv(samplepts1, "C:/Users/Rachel/OneDrive - The University of Montana/Personal/Points/Renne2021_Points.csv")
samplepts1 <- read.csv("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Points/Renne2021_Points.csv", row.names = 1)

# Fix plotseq and etc
head(samplepts1)
samplepts2 <- samplepts1[order(samplepts1$siteID),]

# Now create sequential plot ID's
plotIDs <- rep(c(1:20), nrow(samplepts1)/20)

zeroseq <- rep(c(rep("00",9), rep("0",11)),nrow(samplepts1)/20)

samplepts2$plotseq = plotIDs

samplepts2$PlotID <- paste0(samplepts1$czonename, samplepts1$pzone, "-", zeroseq, samplepts1$plotseq)
row.names(samplepts2) <- seq(1:nrow(samplepts2))

# Save to file:
#write.csv(samplepts2, "C:/Users/Rachel/OneDrive - The University of Montana/Personal/Points/Renne2021_Points.csv")

# Visualize
plot(samplepts1[samplepts1$plotseq == 1, ]$ycoord~samplepts1[samplepts1$plotseq == 1, ]$xcoord, 
     col = c("brown","yellow","purple","grey","red","orange","green","blue","pink")[samplepts1[samplepts1$plotseq == 1, ]$czone], pch = 16, ylim = c(32,49), xlim = c(-125, -100))
map('state', interior = T, add = T)

plot(samplepts1[samplepts1$plotseq < 4, ]$ycoord~samplepts1[samplepts1$plotseq < 4, ]$xcoord, 
     col = c("brown","yellow","purple","grey","red","orange","green","blue","pink")[samplepts1[samplepts1$plotseq < 4, ]$czone], pch = 16, ylim = c(32,49), xlim = c(-125, -100))
map('state', interior = T, add = T)

plot(samplepts1[samplepts1$plotseq < 11, ]$ycoord~samplepts1[samplepts1$plotseq < 11, ]$xcoord, 
     col = c("brown","yellow","purple","grey","red","orange","green","blue","pink")[samplepts1[samplepts1$plotseq < 11, ]$czone], pch = 16, ylim = c(32,49), xlim = c(-125, -100))
map('state', interior = T, add = T)

# Read in climate zones:
czone <- raster("C:/Users/Rachel/OneDrive - The University of Montana/Personal/Sagebrush Sampling Areas/CandidateSamplingAreas_2km_lowslopes_czone.tif")
image(czone, col = c("brown","yellow","purple","grey","red","orange","green","blue","pink"), 
      breaks = c(0,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5), ylim = c(32,49), xlim = c(-125,-100))
map('state', interior = T, add = T)
points(samplepts1[samplepts1$plotseq == 1, ]$ycoord~samplepts1[samplepts1$plotseq == 1, ]$xcoord, 
     col = c("brown","yellow","purple","grey","red","orange","green","blue","pink")[samplepts1[samplepts1$plotseq == 1, ]$czone], pch = 16)
points(samplepts1$ycoord~samplepts1$xcoord, 
       col = c("brown","yellow","purple","grey","red","orange","green","blue","pink")[samplepts1$czone], pch = 16)
points(samplepts1[is.na(samplepts1$czone),]$ycoord[21]~samplepts1[is.na(samplepts1$czone),]$xcoord[21],col = "black", pch = 16)

       