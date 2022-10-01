###########################
# Purpose: Preliminary analyses of plant community
#
#
#
#
#
#
#
###########################

# Bring in data from 2021
plots <- read.csv("C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2021/03_Data_final/plots.csv")
# Keep only relevant columns
plots1 <- plots[,c(1:8)]
names(plots1)[c(1,4,5,6,7,8)] <- c("Plot.Name","Elevation","Date","Auger.depth","Aspect","Slope.Shape")

# Hummocks
humx <- read.csv("C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2021/03_Data_final/hummocks.csv")
# Summarize by plot
humx1 <- data.frame(tapply(humx$mineralsoil, humx$plot, mean))
humx1$Plot.Name <- rownames(humx1)
rownames(humx1) <- 1:nrow(humx1)
names(humx1)[1] <- "avg_soil_cm"
# Change -0.6 value to 0 (no hummock)
humx1[humx1$avg_soil_cm < 0,]$avg_soil_cm <- 0

# Merge onto plots2
plots3 <- merge(plots1, humx1, by = "Plot.Name", all = TRUE)
# Add slope column with NA's
plots3$Percent.slope <- NA

# Par down pall
pall <- read.csv("C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2022/02_Data_final/plots_2022.csv")
pall1 <- pall[,1:8]

# Rearrange plots3 to match pall1
plots4 <- plots3[,c(1:8)]
names(plots4) <- names(pall1)

# Check that names match
table(names(pall1) == names(plots4))

# Combine all plots
allpts <- rbind(plots4, pall1)

##################################
# Get climate data

# Get climate data from daymet & make WL diagrams
library(daymetr)
library(climatol)

# Make day/month list
month <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),
           rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))

# Make dataframes:
ppt = data.frame(plot = allpts$plot, Jan = NA, Feb = NA, Mar = NA, Apr = NA, 
                 May = NA, Jun = NA,Jul = NA, Aug = NA, Sep = NA, Oct = NA, 
                 Nov = NA, Dec = NA)
tmax = ppt
tmin = ppt
tabsmin = ppt

# Run through points 
for (i in 1:nrow(allpts)){
  dmt <- download_daymet(site = "Daymet", lat = allpts$latitude[i], 
                         lon = allpts$longitude[i],start = 1991,end = 2020, 
                         internal = TRUE,silent = TRUE) 
  # Add month column
  dmt$data$month <- month 
  # Create data frame
  ppt[i,2:13] = tapply(dmt$data$prcp..mm.day., dmt$data$month, 
                       FUN = function(x){sum(x)/30})
  tmax[i,2:13] = tapply(dmt$data$tmax..deg.c., dmt$data$month, mean)
  tmin[i,2:13] = tapply(dmt$data$tmin..deg.c., dmt$data$month, mean)
  tabsmin[i,2:13] = tapply(dmt$data$tmin..deg.c., dmt$data$month, min)
}

nrow(allpts)

# Make WL diagrams
for (i in 1:nrow(allpts)){
  png(paste0("C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/WL Diagrams/",allpts$plot[i],".png"),
      width = 5, height = 5, units = "in", res = 300)
  dat <- rbind(ppt[i,], tmax[i,], tmin[i,], tabsmin[i,])
  dat <- dat[,-1]
  diagwl(dat, est = ppt$plot[i], alt = round(allpts$altitude[i]),
         per = "1991-2020")
  dev.off()
}

# Put map and mat onto coords
allpts$map <- apply(ppt[2:13], 1, mean)*10
allpts$mat <- apply(cbind(tmin[,2:13],tmax[,2:13]), 1, mean)

# Save plots to file
#write.csv(allpts, "C:/Users/rache/Dropbox/Doctoral Projects/Projects/Shrub-Grass Relationship Project/Test data merge/allplots.csv",
#          row.names = FALSE)

#allpts <- read.csv("C:/Users/rache/Dropbox/Doctoral Projects/Projects/Shrub-Grass Relationship Project/Test data merge/allplots.csv")

##############################
# Look at data

dir21 <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2021/03_Data_final/"
dir22 <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2022/02_Data_final'

# Make a column for sample year
allpts$year <- c(rep(2021,17),rep(2022,33))

setwd("C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/Preliminary Analysis")

png("SamplePoints_20220912.png", width = 4, height = 4, units = "in", res = 300)
par(mar = c(1,1,1,1))
maps::map("state", regions = c("Wyoming","Nevada","Idaho","Colorado","Arizona",
                               "Utah","Washington","Oregon","Montana","California",
                               "New Mexico"))
points(allpts$latitude~allpts$longitude, col = c("red","blue")[as.factor(allpts$year)], pch = 16)
legend("bottomleft", legend = c("2021","2022"), col = c("red","blue"), pch = 16, bty = "n")
dev.off()

###############################
# Now let's get our one measure per plot for PG, PF, Aherb, Pherb

# Start with 2021 plots
library(dplyr)

# Read in community data
community <- read.csv(file.path(dir21,"community.csv"))

# Read in plant list
plants <- read.csv(file.path(dir21,"codes_edited.csv"))
# Get just codes
plants1 <- plants[,c(1,5:6)]
# Convert Growthform "SE" to "PG"
plants1[!is.na(plants1$GrowthForm) & plants1$GrowthForm == "SE",]$GrowthForm <- "PG"

# Merge codes with community
com <- merge(community, plants1, by.x = "species", by.y = "Symbol", all.x = T)

# Remove observations with no gaps (and corresponding canopy cover microsite)
nogaps <- com[com$Code == "NO GAPS",]

com1 <- com[com$Code != "NO GAPS", ]
for (i in 1:nrow(nogaps)){
  com1 <- com1[-which(com1$plot == nogaps$plot[i] & 
                        com1$transect == nogaps$transect[i] &
                        com1$point == nogaps$point[i]),]
}


# Change from 0,1 to "Open", "Covered"
com1$canopy.code <- NA
com1[com1$canopy == 0,]$canopy.code <- "O"
com1[com1$canopy == 1 | com1$canopy == 0.5,]$canopy.code <- "C"

# Add cover column to com1
com1$cover <- com1$length * com1$width

# Make small cover amounts (0.25 cm2) for rhizomatous grasses
# And make cover 0 when there were no plants observed in that microsite
for (i in 1:nrow(com1)){
  if (com1$number[i] != 0 && is.na(com1$cover[i])){
    com1$cover[i] <- 0.25 
  } 
  if (is.na(com1$GrowthForm[i]) | com1$GrowthForm[i] == "NO PLANTS"){
    com1$cover[i] <- 0
  } 
}


# Get average cover of PG (PG+SE), PF, & AF under and between shrubs:

# Make results DF
results <- data.frame(plot = sort(allpts$plot), PG_microsite_cover_rii = NA, 
                      PG_canopy_cover_rii = NA,
                      PG_microsite_height_rii = NA, PG_canopy_height_rii = NA,
                      PG_microsite_number_rii = NA, PG_canopy_number_rii = NA,
                      PF_microsite_cover_rii = NA, 
                      PF_canopy_cover_rii = NA,
                      PF_microsite_height_rii = NA, PF_canopy_height_rii = NA,
                      PF_microsite_number_rii = NA, PF_canopy_number_rii = NA,
                      Aherb_microsite_cover_rii = NA, 
                      Aherb_canopy_cover_rii = NA,
                      Aherb_microsite_height_rii = NA, Aherb_canopy_height_rii = NA,
                      Aherb_microsite_number_rii = NA, Aherb_canopy_number_rii = NA,
                      Pherb_microsite_cover_rii = NA, 
                      Pherb_canopy_cover_rii = NA,
                      Pherb_microsite_height_rii = NA, Pherb_canopy_height_rii = NA,
                      Pherb_microsite_number_rii = NA, Pherb_canopy_number_rii = NA)

# Pull out individual pfts
PG <- com1[com1$GrowthForm == "PG",]
PF <- com1[com1$GrowthForm == "PF",]
Aherb <- com1[com1$GrowthForm == "AF" | com1$GrowthForm == "AF",]
Pherb <- com1[com1$GrowthForm == "PF"| com1$GrowthForm == "PG",]

pft_summary <- function(thispft){
  # First summarize by microsite
  ms <- thispft %>%
    group_by(plot, transect, point, microsite) %>%
    summarize(num_plants = sum(number, na.rm = T),
              avg_height = mean(height, na.rm = T),
              cover = sum(cover, na.rm = T),
              code = names(table(Code)[which.max(table(Code))]))
  ms <- ms %>%
    group_by(plot, microsite) %>%
    summarize(num_plants = mean(num_plants, na.rm = T),
              avg_height = mean(avg_height, na.rm = T),
              cover = mean(cover, na.rm = T),
              code = names(table(code)[which.max(table(code))]))
  # Add back in microsites that are missing with zeros, etc
  unq <- data.frame(plot = rep(unique(com1$plot),2), microsite = c(rep("C",length(unique(com1$plot))),
                                                                   rep("I",length(unique(com1$plot)))))
  ms <- merge(unq, ms, by = c("plot","microsite"), all = TRUE)
  if (nrow(ms[is.na(ms$num_plants),]) != 0){
    ms[is.na(ms$num_plants),]$num_plants <- 0
  }
  if (nrow(ms[is.na(ms$avg_height),]) != 0){
    ms[is.na(ms$avg_height),]$avg_height <- 0
  }
  if (nrow(ms[is.na(ms$cover),]) != 0){
    ms[is.na(ms$cover),]$cover <- 0
  }
  
  # Then summarize by canopy condition
  ca <- thispft %>%
    group_by(plot, transect, point, canopy.code) %>%
    summarize(num_plants = sum(number, na.rm = T),
              avg_height = mean(height, na.rm = T),
              cover = sum(cover, na.rm = T),
              code = names(table(Code)[which.max(table(Code))]))
  ca <- ca %>%
    group_by(plot, canopy.code) %>%
    summarize(num_plants = mean(num_plants, na.rm = T),
              avg_height = mean(avg_height, na.rm = T),
              cover = mean(cover, na.rm = T),
              code = names(table(code)[which.max(table(code))]))
  unq <- data.frame(plot = rep(unique(com1$plot),2), canopy.code = c(rep("C",length(unique(com1$plot))),
                                                                     rep("O",length(unique(com1$plot)))))
  ca <- merge(unq, ca, by = c("plot","canopy.code"), all = TRUE)
  if (nrow(ca[is.na(ca$num_plants),]) != 0){
    ca[is.na(ca$num_plants),]$num_plants <- 0
  }
  if (nrow(ca[is.na(ca$avg_height),]) != 0){
    ca[is.na(ca$avg_height),]$avg_height <- 0
  }
  if (nrow(ca[is.na(ca$cover),]) != 0){
    ca[is.na(ca$cover),]$cover <- 0
  }
  
  return(list(ms,ca))
}

PG_summary <- pft_summary(PG)
PF_summary <- pft_summary(PF)
Aherb_summary <- pft_summary(Aherb)
Pherb_summary <- pft_summary(Pherb)

##########################################################################
# Read back in hand-entered codes
codes <- read.csv(file.path(dir22, "codes_2022.csv"))

# Read in pcomm from 2022
pcomm1 <- read.csv(file.path(dir22,"community_2022.csv"))

# Merge codes with pcomm
pcomm2 <- merge(pcomm1, codes[,c(1,5,6)], by.x = "species", by.y = "code", all.x = TRUE)



# In Number column, set to 1 if NA
pcomm2[is.na(pcomm2$number) & pcomm2$species != "NO PLANTS", ]$number <- 1
# Set Number to 0 where "species" is "NO PLANTS"
pcomm2[pcomm2$species == "NO PLANTS",]$number <- 0
# Calculate Number of plants
sum(na.omit(pcomm2$number))
# 11068

# Remove observations with no gaps (and corresponding canopy cover microsite)
nogaps <- pcomm2[pcomm2$species == "NO GAPS >60 cm",]

pcomm3 <- pcomm2[pcomm2$species != "NO GAPS >60 cm", ]
for (i in 1:nrow(nogaps)){
  pcomm3 <- pcomm3[-which(pcomm3$plot == nogaps$plot[i] & 
                            pcomm3$transect == nogaps$transect[i] &
                            pcomm3$point == nogaps$point[i]),]
}

# Make best guess about missing Canopy.condition
pcomm3[pcomm3$plot == "GB6-011" & pcomm3$transect == 2 & pcomm3$point == 5 & pcomm3$microsite == "canopy",]
pcomm3[is.na(pcomm3$canopy),]$canopy <- 1


# Change from 0,1 to "Open", "Covered"
pcomm3$canopy.code <- "NONE"
pcomm3[pcomm3$canopy == 0,]$canopy.code <- "O"
pcomm3[pcomm3$canopy == 1 | pcomm3$canopy == 0.5,]$canopy.code <- "C"

# Add cover column to pcomm3
pcomm3$cover <- pcomm3$length * pcomm3$width

# Make small cover amounts (0.25 cm2) for rhizomatous grasses
# And make cover 0 when there were no plants observed in that microsite
for (i in 1:nrow(pcomm3)){
  if (pcomm3$number[i] != 0 && is.na(pcomm3$cover[i])){
    pcomm3$cover[i] <- 0.25 
  } 
  if (is.na(pcomm3$GrowthForm[i]) | pcomm3$GrowthForm[i] == "NO PLANTS"){
    pcomm3$cover[i] <- 0
  } 
}

# Change microsite coding
pcomm3$microsite <- as.character(pcomm3$microsite)
pcomm3$microsite <- gsub("interspace","I",pcomm3$microsite)
pcomm3$microsite <- gsub("canopy","C",pcomm3$microsite)

# Fix "NO PLANTS"
pcomm3[pcomm3$species == "NO PLANTS",]$Code <- "NO PLANTS"

# Pull out individual pfts
PG1 <- pcomm3[!is.na(pcomm3$GrowthForm) & pcomm3$GrowthForm == "PG",]
PF1 <- pcomm3[!is.na(pcomm3$GrowthForm) & pcomm3$GrowthForm == "PF",]
Aherb1 <- pcomm3[!is.na(pcomm3$GrowthForm) & (pcomm3$GrowthForm == "AF" | pcomm3$GrowthForm == "AG"),]
Pherb1 <- pcomm3[!is.na(pcomm3$GrowthForm) & (pcomm3$GrowthForm == "PF"| pcomm3$GrowthForm == "PG"),]

# Troubleshooting example
thispft <- Aherb1

pft_summary1 <- function(thispft){
  # First summarize by microsite
  ms <- thispft %>%
    group_by(plot, transect, point, microsite) %>%
    summarize(num_plants = sum(number, na.rm = T),
              avg_height = mean(height, na.rm = T),
              cover = sum(cover, na.rm = T),
              code = NA)#names(table(Code)[which.max(table(Code))]))
  ms <- ms %>%
    group_by(plot, microsite) %>%
    summarize(num_plants = mean(num_plants, na.rm = T),
              avg_height = mean(avg_height, na.rm = T),
              cover = mean(cover, na.rm = T),
              code = NA)#names(table(code)[which.max(table(code))]))
  # Add back in microsites that are missing with zeros, etc
  unq <- data.frame(plot = rep(unique(pcomm3$plot),2), microsite = c(rep("C",length(unique(pcomm3$plot))),
                                                                     rep("I",length(unique(pcomm3$plot)))))
  ms <- merge(unq, ms, by = c("plot","microsite"), all = TRUE)
  if (nrow(ms[is.na(ms$num_plants),]) != 0){
    ms[is.na(ms$num_plants),]$num_plants <- 0
  }
  if (nrow(ms[is.na(ms$avg_height),]) != 0){
    ms[is.na(ms$avg_height),]$avg_height <- 0
  }
  if (nrow(ms[is.na(ms$cover),]) != 0){
    ms[is.na(ms$cover),]$cover <- 0
  }
  
  # Then summarize by canopy condition
  ca <- thispft %>%
    group_by(plot, transect, point, canopy.code) %>%
    summarize(num_plants = sum(number, na.rm = T),
              avg_height = mean(height, na.rm = T),
              cover = sum(cover, na.rm = T),
              code = NA)#names(table(Code)[which.max(table(Code))]))
  ca <- ca %>%
    group_by(plot, canopy.code) %>%
    summarize(num_plants = mean(num_plants, na.rm = T),
              avg_height = mean(avg_height, na.rm = T),
              cover = mean(cover, na.rm = T),
              code = NA)#names(table(code)[which.max(table(code))]))
  unq <- data.frame(plot = rep(unique(pcomm3$plot),2), canopy.code = c(rep("C",length(unique(pcomm3$plot))),
                                                                       rep("O",length(unique(pcomm3$plot)))))
  ca <- merge(unq, ca, by = c("plot","canopy.code"), all = TRUE)
  if (nrow(ca[is.na(ca$num_plants),]) != 0){
    ca[is.na(ca$num_plants),]$num_plants <- 0
  }
  if (nrow(ca[is.na(ca$avg_height),]) != 0){
    ca[is.na(ca$avg_height),]$avg_height <- 0
  }
  if (nrow(ca[is.na(ca$cover),]) != 0){
    ca[is.na(ca$cover),]$cover <- 0
  }
  
  return(list(ms,ca))
}

PG_summary1 <- pft_summary1(PG1)
PF_summary1 <- pft_summary1(PF1)
Aherb_summary1 <- pft_summary1(Aherb1)
Pherb_summary1 <- pft_summary1(Pherb1)

##################################################
# Combine all plots (make function)
combine_plots <- function(df1, df2){
  df1a <- cbind(df1[[1]], df1[[2]][,-1])
  names(df1a) <- c("Plot","Microsite","ms_num_plants","ms_avg_height","ms_cover","ms_code",
                   "canopy.code","ca_num_plants","ca_avg_height","ca_cover","ca_code" )
  df2a <- cbind(df2[[1]], df2[[2]][,-1])
  names(df2a) <- c("Plot","Microsite","ms_num_plants","ms_avg_height","ms_cover","ms_code",
                   "canopy.code","ca_num_plants","ca_avg_height","ca_cover","ca_code" )
  return(rbind(df1a, df2a))
}

# Combine all plots with function
PGall <- combine_plots(PG_summary, PG_summary1)
PFall <- combine_plots(PF_summary, PF_summary1)
Pherball <- combine_plots(Pherb_summary, Pherb_summary1)
Aherball <- combine_plots(Aherb_summary, Aherb_summary1)

# Add RII's into results
for (i in 1:nrow(results)){
  # PG
  x <- PGall[PGall$Plot == results$plot[i],]
  results$PG_microsite_cover_rii[i] <- (x$ms_cover[which(x$Microsite == "C")] - x$ms_cover[which(x$Microsite == "I")])/
    (x$ms_cover[which(x$Microsite == "C")] + x$ms_cover[which(x$Microsite == "I")])
  results$PG_canopy_cover_rii[i] <- (x$ms_cover[which(x$canopy.code == "C")] - x$ms_cover[which(x$canopy.code == "O")])/
    (x$ms_cover[which(x$canopy.code == "C")] + x$ms_cover[which(x$canopy.code == "O")])
  results$PG_microsite_height_rii[i] <- (x$ms_avg_height[which(x$Microsite == "C")] - x$ms_avg_height[which(x$Microsite == "I")])/
    (x$ms_avg_height[which(x$Microsite == "C")] + x$ms_avg_height[which(x$Microsite == "I")])
  results$PG_canopy_height_rii[i] <- (x$ca_avg_height[which(x$canopy.code == "C")] - x$ca_avg_height[which(x$canopy.code == "O")])/
    (x$ca_avg_height[which(x$canopy.code == "C")] + x$ca_avg_height[which(x$canopy.code == "O")])
  results$PG_microsite_number_rii[i] <- (x$ms_num_plants[which(x$Microsite == "C")] - x$ms_num_plants[which(x$Microsite == "I")])/
    (x$ms_num_plants[which(x$Microsite == "C")] + x$ms_num_plants[which(x$Microsite == "I")])
  results$PG_canopy_number_rii[i] <- (x$ms_num_plants[which(x$canopy.code == "C")] - x$ms_num_plants[which(x$canopy.code == "O")])/
    (x$ms_num_plants[which(x$canopy.code == "C")] + x$ms_num_plants[which(x$canopy.code == "O")])
  # PF
  x <- PFall[PFall$Plot == results$plot[i],]
  results$PF_microsite_cover_rii[i] <- (x$ms_cover[which(x$Microsite == "C")] - x$ms_cover[which(x$Microsite == "I")])/
    (x$ms_cover[which(x$Microsite == "C")] + x$ms_cover[which(x$Microsite == "I")])
  results$PF_canopy_cover_rii[i] <- (x$ms_cover[which(x$canopy.code == "C")] - x$ms_cover[which(x$canopy.code == "O")])/
    (x$ms_cover[which(x$canopy.code == "C")] + x$ms_cover[which(x$canopy.code == "O")])
  results$PF_microsite_height_rii[i] <- (x$ms_avg_height[which(x$Microsite == "C")] - x$ms_avg_height[which(x$Microsite == "I")])/
    (x$ms_avg_height[which(x$Microsite == "C")] + x$ms_avg_height[which(x$Microsite == "I")])
  results$PF_canopy_height_rii[i] <- (x$ca_avg_height[which(x$canopy.code == "C")] - x$ca_avg_height[which(x$canopy.code == "O")])/
    (x$ca_avg_height[which(x$canopy.code == "C")] + x$ca_avg_height[which(x$canopy.code == "O")])
  results$PF_microsite_number_rii[i] <- (x$ms_num_plants[which(x$Microsite == "C")] - x$ms_num_plants[which(x$Microsite == "I")])/
    (x$ms_num_plants[which(x$Microsite == "C")] + x$ms_num_plants[which(x$Microsite == "I")])
  results$PF_canopy_number_rii[i] <- (x$ms_num_plants[which(x$canopy.code == "C")] - x$ms_num_plants[which(x$canopy.code == "O")])/
    (x$ms_num_plants[which(x$canopy.code == "C")] + x$ms_num_plants[which(x$canopy.code == "O")])
  # Pherb
  x <- Pherball[Pherball$Plot == results$plot[i],]
  results$Pherb_microsite_cover_rii[i] <- (x$ms_cover[which(x$Microsite == "C")] - x$ms_cover[which(x$Microsite == "I")])/
    (x$ms_cover[which(x$Microsite == "C")] + x$ms_cover[which(x$Microsite == "I")])
  results$Pherb_canopy_cover_rii[i] <- (x$ms_cover[which(x$canopy.code == "C")] - x$ms_cover[which(x$canopy.code == "O")])/
    (x$ms_cover[which(x$canopy.code == "C")] + x$ms_cover[which(x$canopy.code == "O")])
  results$Pherb_microsite_height_rii[i] <- (x$ms_avg_height[which(x$Microsite == "C")] - x$ms_avg_height[which(x$Microsite == "I")])/
    (x$ms_avg_height[which(x$Microsite == "C")] + x$ms_avg_height[which(x$Microsite == "I")])
  results$Pherb_canopy_height_rii[i] <- (x$ca_avg_height[which(x$canopy.code == "C")] - x$ca_avg_height[which(x$canopy.code == "O")])/
    (x$ca_avg_height[which(x$canopy.code == "C")] + x$ca_avg_height[which(x$canopy.code == "O")])
  results$Pherb_microsite_number_rii[i] <- (x$ms_num_plants[which(x$Microsite == "C")] - x$ms_num_plants[which(x$Microsite == "I")])/
    (x$ms_num_plants[which(x$Microsite == "C")] + x$ms_num_plants[which(x$Microsite == "I")])
  results$Pherb_canopy_number_rii[i] <- (x$ms_num_plants[which(x$canopy.code == "C")] - x$ms_num_plants[which(x$canopy.code == "O")])/
    (x$ms_num_plants[which(x$canopy.code == "C")] + x$ms_num_plants[which(x$canopy.code == "O")])
  # Aherb
  x <- Aherball[Aherball$Plot == results$plot[i],]
  results$Aherb_microsite_cover_rii[i] <- (x$ms_cover[which(x$Microsite == "C")] - x$ms_cover[which(x$Microsite == "I")])/
    (x$ms_cover[which(x$Microsite == "C")] + x$ms_cover[which(x$Microsite == "I")])
  results$Aherb_canopy_cover_rii[i] <- (x$ms_cover[which(x$canopy.code == "C")] - x$ms_cover[which(x$canopy.code == "O")])/
    (x$ms_cover[which(x$canopy.code == "C")] + x$ms_cover[which(x$canopy.code == "O")])
  results$Aherb_microsite_height_rii[i] <- (x$ms_avg_height[which(x$Microsite == "C")] - x$ms_avg_height[which(x$Microsite == "I")])/
    (x$ms_avg_height[which(x$Microsite == "C")] + x$ms_avg_height[which(x$Microsite == "I")])
  results$Aherb_canopy_height_rii[i] <- (x$ca_avg_height[which(x$canopy.code == "C")] - x$ca_avg_height[which(x$canopy.code == "O")])/
    (x$ca_avg_height[which(x$canopy.code == "C")] + x$ca_avg_height[which(x$canopy.code == "O")])
  results$Aherb_microsite_number_rii[i] <- (x$ms_num_plants[which(x$Microsite == "C")] - x$ms_num_plants[which(x$Microsite == "I")])/
    (x$ms_num_plants[which(x$Microsite == "C")] + x$ms_num_plants[which(x$Microsite == "I")])
  results$Aherb_canopy_number_rii[i] <- (x$ms_num_plants[which(x$canopy.code == "C")] - x$ms_num_plants[which(x$canopy.code == "O")])/
    (x$ms_num_plants[which(x$canopy.code == "C")] + x$ms_num_plants[which(x$canopy.code == "O")])
}

# Merge results with allpts
allpts1 <- merge(allpts, results, by = "plot", all = TRUE)

# Add in a column for climate zone
allpts1$czone <- gsub("GP.*","GP", allpts1$plot)
allpts1$czone <- gsub("GB.*","GB", allpts1$czone)
allpts1$czone <- gsub("DryGB.*","DryGB", allpts1$czone)
allpts1$czone <- gsub("CP.*","CP", allpts1$czone)
allpts1$czone <- gsub("PNWMtn.*","PNWMtn", allpts1$czone)
allpts1$czone <- gsub("RkyMtn.*","RkyMtn", allpts1$czone)
allpts1$czone <- gsub("WYBas.*","WYBas", allpts1$czone)

# Make color ramp
zonecols <- colorRampPalette(c('orange','yellow','blue','purple','brown','pink','forestgreen'), alpha = TRUE)

# png("Hummocks_vs_map_20220719.png", width = 4, height = 4, units = "in", res = 300)
# par(mar= c(2,2,2,1), mgp = c(1,0.1,0), tcl =0.1)
# plot(allpts1$avg_soil_cm~allpts1$map,
#      col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
#      xlab = "Mean annual precipitation (mm)",
#      ylab = "Average hummock depth (cm)",
#      main = "Hummocks")
# dev.off()

# Look at climate variables
png("Mat_vs_map_20220912.png", width = 4, height = 4, units = "in", res = 300)
par(mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$map~allpts1$mat, main = "Climate of plots",
     xlab = "Mean annual temperature (C)",
     ylab = "Mean annual precipitation (mm)",
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16)
dev.off()

#######################################
# PG
png("PG_cover_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$PG_microsite_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$PG_microsite_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$PG_canopy_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$PG_canopy_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()


png("PG_height_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$PG_microsite_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$PG_microsite_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$PG_canopy_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$PG_canopy_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

png("PG_number_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$PG_microsite_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$PG_microsite_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$PG_canopy_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$PG_canopy_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

#######################################
# PF
png("PF_cover_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$PF_microsite_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$PF_microsite_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$PF_canopy_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$PF_canopy_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()


png("PF_height_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$PF_microsite_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$PF_microsite_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$PF_canopy_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$PF_canopy_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

png("PF_number_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$PF_microsite_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$PF_microsite_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$PF_canopy_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$PF_canopy_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

#######################################
# Pherb
png("Pherb_cover_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$Pherb_microsite_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$Pherb_microsite_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$Pherb_canopy_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$Pherb_canopy_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()


png("Pherb_height_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$Pherb_microsite_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$Pherb_microsite_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$Pherb_canopy_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$Pherb_canopy_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

png("Pherb_number_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$Pherb_microsite_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$Pherb_microsite_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$Pherb_canopy_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$Pherb_canopy_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

#######################################
# Aherb
png("Aherb_cover_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$Aherb_microsite_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$Aherb_microsite_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$Aherb_canopy_cover_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$Aherb_canopy_cover_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()


png("Aherb_height_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$Aherb_microsite_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$Aherb_microsite_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$Aherb_canopy_height_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$Aherb_canopy_height_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()

png("Aherb_number_RIIs.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2,2), mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(allpts1$Aherb_microsite_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Microsites")
abline(h = 0, lty = 3)
plot(allpts1$Aherb_microsite_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Microsites")
abline(h = 0, lty = 3)

plot(allpts1$Aherb_canopy_number_rii~allpts1$map,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual precipitation (mm)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
plot(allpts1$Aherb_canopy_number_rii~allpts1$mat,
     col = zonecols(7)[as.factor(allpts1$czone)], pch = 16,
     ylim = c(-1,1), ylab = "RII (Competitive to Facultative)",
     xlab = "Mean annual temperature (C)",
     main = "Canopy condition")
abline(h = 0, lty = 3)
dev.off()


