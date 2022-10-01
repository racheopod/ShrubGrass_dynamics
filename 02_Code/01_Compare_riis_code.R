##################################
# Purpose: Calculate and compare relative interaction index (rii) 
#          between community and lpi measurements
#
#
# Rachel R. Renne
# September 21, 2022
#################################

# Set up directories
datadir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data'
codedir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/'

# Read in code files with functions to calculate RIIs
source(file.path(codedir, "lpi.R"))
source(file.path(codedir, "plantcommunity.R"))

# Read in plots file
plots <- read.csv(file.path(datadir, "plots.csv"))

# Remove notes column
plots1 <- plots[,-13]

# Add a climate zone column
plots1$czone <- plots$plot
czones <- list(c("CP.*","CP"),c("DryGB.*","DryGB"),c("GB.*","GB"),
            c("GP.*","GP"),c("PNWMtn.*","PNWMtn"),c("WYBas.*","WYBas"),
            c("RkyMtn.*","RkyMtn"))
for (i in 1:length(czones)){
  plots1$czone <- gsub(czones[[i]][1],czones[[i]][2],plots1$czone)
}

# Make vector of colors to correspond to czones
zonecols <- colorRampPalette(c('orange','yellow','blue','purple','brown','pink','forestgreen'), alpha = TRUE)

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
cover <- lpicover(lpi, n_layers = 5, separate_dead = FALSE,
                   nonplants = c("VL","HL","WL"),
                   soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"),
                   pres.abs_abund = FALSE)

# Pull out just sagebrush cover
sage <- cover[cover$species %in% c("ARTRW8","ARTRT","ARTRV","ARNO4","ARCA13","ARARA"),]

# Change name of cover column to reflect sagebrush cover
names(sage)[2:3] <- c("sagebrush_species","sagebrush_cover")

# Merge with plots
plots2 <- merge(plots1, sage, by = "plot")

# Calculate riis of cover based on microsites (under vs between shrubs)
lpi_riis <- lpi_microsite(lpi, ftypes, n_layers = 5,
                      soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"))

# Change rii name for clarity during comparison 
names(lpi_riis)[12] <- "lpi_cover_riis"

# Remove cases where interspace and canopy cover both equal 0
lpi_riis1 <- lpi_riis[lpi_riis$canopy_count > 0 & lpi_riis$interspace_count > 0 & lpi_riis$species != "N", ]

# Merge lpi_riis with plots
lpi_riis2 <- merge(plots2, lpi_riis1[,c(1,2,12)])

# Read in community data
pcomm <- read.csv(file.path(datadir, "community.csv"))

# Calculate both microsite and canopy condition values
ms_and_ca_calcs <- pcomm_microsite(pcomm, ftype, can_edge = FALSE)

# Separate out each one
microsites <- ms_and_ca_calcs[[1]]
canopy_condition <- ms_and_ca_calcs[[2]]

# Remove microsite records where there are 0 in canopy and microsite interspaces
microsites1 <- microsites[microsites$canopy_n_ind > 0 & microsites$interspace_n_ind > 0, ]

# Calculate some RIIs
microsites1$pcomm_cover_rii <- calculate_rii(canopy=microsites1$canopy_avg_cover, 
                                            interspace=microsites1$interspace_avg_cover)

microsites1$pcomm_number_rii <- calculate_rii(canopy=microsites1$canopy_avg_number, 
                                            interspace=microsites1$interspace_avg_number)

microsites1$pcomm_size_rii <- calculate_rii(canopy=microsites1$canopy_size,
                                           interspace=microsites1$interspace_size)

microsites1$pcomm_height_rii <- calculate_rii(canopy=microsites1$canopy_height,
                                             interspace=microsites1$interspace_height)
names(microsites1)[2] <- "species"

# Merge with lpi_riis
riis <- merge(lpi_riis1[,c(1,2,12)], microsites1[,c(1,2,5:8,11:14,17:20)], by = c("plot","species"), all = TRUE)

# Merge with plots
riis1 <- merge(plots2, riis, by = "plot")

# Pull out individual ftypes for comparison
shrubs <- riis1[riis1$species == "SH",]
pgrass <- riis1[riis1$species == "PG",]
agrass <- riis1[riis1$species == "AG",]
pforb <- riis1[riis1$species == "PF",]
aherb <- riis1[riis1$species %in% c("AG","AF"),]
pherb <- riis1[riis1$species %in% c("PG","PF"),]
allherb <- riis1[riis1$species %in% c("AG","AF","PG","PF"),]

# Summarize over plots for each
combineFT <- function(df){
  `%>%` <- magrittr::`%>%`
  df <- df %>% 
    dplyr::group_by(plot) %>%
    dplyr::summarize(canopy_avg_number = sum(canopy_avg_number, na.rm = T),
                     interspace_avg_number = sum(interspace_avg_number, na.rm = T),
                     canopy_avg_cover = sum(canopy_avg_cover, na.rm = T),      
                     interspace_avg_cover = sum(interspace_avg_cover, na.rm = T),  
                     canopy_size = sum(canopy_size, na.rm = T),          
                     interspace_size = sum(interspace_size, na.rm = T),       
                     canopy_height = sum(canopy_height, na.rm = T),        
                     interspace_height = sum(interspace_height, na.rm = T))
  # Calculate some RIIs
  df$pcomm_cover_rii <- calculate_rii(canopy=df$canopy_avg_cover, 
                                               interspace=df$interspace_avg_cover)
  
  df$pcomm_number_rii <- calculate_rii(canopy=df$canopy_avg_number, 
                                                interspace=df$interspace_avg_number)
  
  df$pcomm_size_rii <- calculate_rii(canopy=df$canopy_size,
                                              interspace=df$interspace_size)
  
  df$pcomm_height_rii <- calculate_rii(canopy=df$canopy_height,
                                                interspace=df$interspace_height)
  df <- merge(plots2, df, by = "plot")
  return(df)
} 

# Combine fts
aherb <- combineFT(aherb)
pherb <- combineFT(pherb)
allherb <- combineFT(allherb)

# Start with pgrass

compare_riis <- function(rii2, rii1, czones, rii2_name, rii1_name, ftype){
  par(mar=c(2,2,2,1),mgp=c(1,0.1,0),tcl=0.1)
  plot(rii2~rii1, main = paste0("Comparing RIIs: ", ftype),
       xlab = rii1_name, ylab = rii2_name,
       col = zonecols(7)[as.factor(czones)], pch = 16,
       xlim = c(-1,1), ylim = c(-1,1))
  abline(h=0,lty=3)
  abline(v=0,lty=3)
}

compare_riis(pgrass$pcomm_cover_rii, pgrass$lpi_cover_riis, pgrass$czone,
             "Plant Community average cover", "LPI Cover","PGrass")
compare_riis(pgrass$pcomm_size_rii, pgrass$lpi_cover_riis, pgrass$czone,
             "Plant Community average size", "LPI Cover", "PGrass")
compare_riis(pgrass$pcomm_height_rii, pgrass$lpi_cover_riis, pgrass$czone,
             "Plant Community average height", "LPI Cover", "PGrass")
compare_riis(pgrass$pcomm_number_rii, pgrass$lpi_cover_riis, pgrass$czone,
             "Plant Community average number/600cm2", "LPI Cover", "PGrass")

par(mar=c(2,2,2,1),mgp=c(1,0.1,0),tcl=0.1)
plot(pgrass$lpi_cover_riis~pgrass$map, ylim = c(-1,1),
     col = zonecols(7)[as.factor(pgrass$czone)], pch = 16)
abline(h=0,lty=3)
plot(pgrass$pcomm_cover_rii~pgrass$map, ylim = c(-1,1),
     col = zonecols(7)[as.factor(pgrass$czone)], pch = 16)
abline(h=0,lty=3)


# Look at PForb
compare_riis(pforb$pcomm_cover_rii, pforb$lpi_cover_riis, pforb$czone,
             "Plant Community average cover", "LPI Cover","PForb")
compare_riis(pforb$pcomm_size_rii, pforb$lpi_cover_riis, pforb$czone,
             "Plant Community average size", "LPI Cover", "PForb")
compare_riis(pforb$pcomm_height_rii, pforb$lpi_cover_riis, pforb$czone,
             "Plant Community average height", "LPI Cover", "PForb")
compare_riis(pforb$pcomm_number_rii, pforb$lpi_cover_riis, pforb$czone,
             "Plant Community average number/600cm2", "LPI Cover", "PForb")



####################Make some RII plots:

riiplot <- function(x,y, xlab,main, col, cex,add.plot = TRUE, plotrow = length(y),...){
  par(mar=c(2,2,2,1),mgp=c(1,0.1,0),tcl=0.1)
  plot(y~x, pch = 16, ylim = c(-1,1),
       col = col,
       cex = cex,
       main = main,
       xlab = xlab,
       ylab = "RII (competitive to facultative)")
  if (add.plot){
    points(x = x[plotrow], y = y[plotrow],
           col = "red", pch = 16, cex = 0.7)
  }
  abline(h=0,lty=3)
}

fig_path <- "C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/05_Data_Analysis/03_Figures/preliminary"

png(file.path(fig_path,"Aherb.png"), width = 7, height = 3, units = "in", res = 300)
par(mfrow = c(1,2))
riiplot(y = aherb$pcomm_number_rii, x = aherb$map, xlab = "Mean annual precipitation (mm)",
        main = "Aherb: average number", cex = 0.7, col = zonecols(7)[as.factor(aherb$czone)])
riiplot(y = aherb$pcomm_height_rii, x = aherb$map, xlab = "Mean annual precipitation (mm)",
        main = "Aherb: average height", cex = 0.7, col = zonecols(7)[as.factor(aherb$czone)])
dev.off()

png(file.path(fig_path,"Agrass.png"), width = 7, height = 3, units = "in", res = 300)
par(mfrow = c(1,2))
riiplot(y = agrass$pcomm_number_rii, x = agrass$map, xlab = "Mean annual precipitation (mm)",
        main = "Agrass: average number", cex = 0.7, col = zonecols(7)[as.factor(agrass$czone)],
        add.plot = FALSE)
riiplot(y = agrass$pcomm_height_rii, x = agrass$map, xlab = "Mean annual precipitation (mm)",
        main = "Agrass: average height", cex = 0.7, col = zonecols(7)[as.factor(agrass$czone)],
        add.plot = FALSE)
dev.off()

png(file.path(fig_path,"Pherb.png"), width = 7, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
riiplot(y = pherb$pcomm_number_rii, x = pherb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pherb: average number", cex = 0.7, col = zonecols(7)[as.factor(pherb$czone)])
riiplot(y = pherb$pcomm_height_rii, x = pherb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pherb: average height", cex = 0.7, col = zonecols(7)[as.factor(pherb$czone)])
riiplot(y = pherb$pcomm_cover_rii, x = pherb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pherb: cover", cex = 0.7, col = zonecols(7)[as.factor(pherb$czone)])
riiplot(y = pherb$pcomm_size_rii, x = pherb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pherb: average size", cex = 0.7, col = zonecols(7)[as.factor(pherb$czone)])
dev.off()

png(file.path(fig_path,"Pgrass.png"), width = 7, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
riiplot(y = pgrass$pcomm_number_rii, x = pgrass$map, xlab = "Mean annual precipitation (mm)",
        main = "Pgrass: average number", cex = 0.7, col = zonecols(7)[as.factor(pgrass$czone)])
riiplot(y = pgrass$pcomm_height_rii, x = pgrass$map, xlab = "Mean annual precipitation (mm)",
        main = "Pgrass: average height", cex = 0.7, col = zonecols(7)[as.factor(pgrass$czone)])
riiplot(y = pgrass$pcomm_cover_rii, x = pgrass$map, xlab = "Mean annual precipitation (mm)",
        main = "Pgrass: cover", cex = 0.7, col = zonecols(7)[as.factor(pgrass$czone)])
riiplot(y = pgrass$pcomm_size_rii, x = pgrass$map, xlab = "Mean annual precipitation (mm)",
        main = "Pgrass: average size", cex = 0.7, col = zonecols(7)[as.factor(pgrass$czone)])
dev.off()

png(file.path(fig_path,"Pforb.png"), width = 7, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
riiplot(y = pforb$pcomm_number_rii, x = pforb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pforb: average number", cex = 0.7, col = zonecols(7)[as.factor(pforb$czone)])
riiplot(y = pforb$pcomm_height_rii, x = pforb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pforb: average height", cex = 0.7, col = zonecols(7)[as.factor(pforb$czone)])
riiplot(y = pforb$pcomm_cover_rii, x = pforb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pforb: cover", cex = 0.7, col = zonecols(7)[as.factor(pforb$czone)])
riiplot(y = pforb$pcomm_size_rii, x = pforb$map, xlab = "Mean annual precipitation (mm)",
        main = "Pforb: average size", cex = 0.7, col = zonecols(7)[as.factor(pforb$czone)])
dev.off()

png(file.path(fig_path,"Allherb.png"), width = 7, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
riiplot(y = allherb$pcomm_number_rii, x = allherb$map, xlab = "Mean annual precipitation (mm)",
        main = "AllHerb: average number", cex = 0.7, col = zonecols(7)[as.factor(allherb$czone)])
riiplot(y = allherb$pcomm_height_rii, x = allherb$map, xlab = "Mean annual precipitation (mm)",
        main = "AllHerb: average height", cex = 0.7, col = zonecols(7)[as.factor(allherb$czone)])
riiplot(y = allherb$pcomm_cover_rii, x = allherb$map, xlab = "Mean annual precipitation (mm)",
        main = "AllHerb: cover", cex = 0.7, col = zonecols(7)[as.factor(allherb$czone)])
riiplot(y = allherb$pcomm_size_rii, x = allherb$map, xlab = "Mean annual precipitation (mm)",
        main = "AllHerb: average size", cex = 0.7, col = zonecols(7)[as.factor(allherb$czone)])
dev.off()

png(file.path(fig_path, "Climatespace.png"), width = 3.5, height = 3, units = "in",
    res = 200)
par(mar=c(2,2,1,1),mgp=c(1,0.1,0),tcl=0.1)
plot(x = plots2$map, y = plots2$mat, col = zonecols(7)[as.factor(plots2$czone)],
     pch = 16, xlab = "Mean annual precipitation (mm)", cex = 0.7,
     ylab = "Mean annual temperature")
points(x = plots[plots$plot == "WYBasX-001",]$map, plots[plots$plot == "WYBasX-001",]$mat,
       cex = 0.7, pch = 16, col = "red")
dev.off()
