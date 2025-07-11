#' lpicover_ftypecorrected: Extract plant functional type cover from line-point 
#' intercept data, where functional types are combined at the point level.
#' 
#' 
#' @param lpi A dataframe of line-point intercept observations. This dataframe
#'  should include columns for `plot`, `transect`, and `point`, followed by  
#'  columns named `toplayer`,`layer1`,`layer2`,`layer3`,... up to `n_layers`-1. 
#'  If `separate_dead` is `TRUE`, there should be additional columns corresponding 
#'  to each layer column (named following the pattern `toplayer_dead`,`layer1_dead`,
#'  `layer2_dead`, etc.) that designate `1` if the plant part in that layer was 
#'  recorded as dead, and `0` otherwise.
#' @param n_layers Number of layers.
#' @param soilsurface vector. A character vector of all non-plant soil surface 
#'  codes present in the data (default = c(`BR`,`BY`,`CB`,`EL`,`GR`,`LC`,`M`,`S`,`ST`).
#' @param pres.abs_abund boolean. `TRUE` indicates a output will be a dataframe
#'  with rows corresponding to plots and columns corresponding to all species in 
#'  the dataset so that species not detected at a given plot in lpi will be 
#'  assigned a `0` (while fractional cover will be given for species present in 
#'  lpi). `FALSE` indicates that output will be a dataframe with three
#'  columns: `plot`,`species` (i.e. cover type), and `cover` (fractional). default
#'  is `FALSE`.
#' @return data frame of cover values for each plot, with three columns: `plot`,
#'  `species` (i.e. cover type), and `cover` (fractional). (Created 20250114).
#'  
#' @details This deals with the issue of combining individual species cover into 
#'  functional type cover. For example, if a pin hits two POSE plants, we only
#'  record POSE one time, but if it hits POSE and ELEL5, we record both. When
#'  summing by functional type, the former case gets one PG hit but the latter case
#'  gets two. This function combines functional types per point, so both cases would
#'  get a single hit. This produces a cover estimate that is easier to interpret
#'  than the common approach of combining species cover by functional type.
#'  
#' @examples
#' # File path to final, cleaned data
#' cleandir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data'
#' 
#' # Read in lpi
#' lpi <- read.csv(file.path(cleandir,"lpi.csv"))
#' 
#' # Create the ftype data frame
#' codes <- read.csv(file.path(cleandir,"codes.csv"))
#' codes1 <- codes[,c("code","GrowthForm")]
#' names(codes1) <- c("species","ftype")
#' # Remove "S" (soil) rows
#' codes2 <- codes1[!(codes1$ftype == "S"),]
#' ftype = codes2

#' # Calculate cover
#' cover <- lpicover_ftypecorrected(lpi, ftypes, n_layers = 5,
#'                                 soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"))

#' par(mar = c(2,2,1,1), mgp = c(1,0.1,0), tcl = 0.1)
#' plot(cover[cover$species == "PF",]$cover ~ 
#'        cover[cover$plot %in% cover[cover$species == "PF",]$plot & cover$species == "SH",]$cover,
#'      ylab = "Perennial forb cover",
#'      xlab = "Shrub cover", xlim = c(0,1), ylim = c(0,1))


lpicover_ftypecorrected <- function(lpi, ftype, n_layers = 5, shrub_code = "SH",
                                    soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST")){
  # Check for soil surface plants that aren't found in layers above
  basalcheck <- lpi[!(lpi$soilsurface %in% soilsurface),]
  layers <- c("toplayer",paste0("layer",1:(n_layers-1)))
  if (sum(apply(basalcheck,1,function(x){!(x["soilsurface"] %in% x[layers])})) > 0){
    warning("One or more plant codes in soilsurface layer may not be present in upper layers.")
  } else {rm(basalcheck)} # Remove basalcheck if no problem
  
  # Verify that there are no values == "" (empty) in top layer or soil surface
  if (sum(lpi$toplayer == "") > 0){
    stop(paste0("'toplayer' is missing for ",sum(lpi$toplayer == "")," records."))
  }
  if (sum(lpi$soilsurface == "") > 0){
    stop(paste0("'soilsurface' is missing for ",sum(lpi$soilsurface == "")," records."))
  }
  
  # Create points/plot named vector
  pts_plt <- data.frame(table(lpi$plot))
  names(pts_plt) <- c("plot","points")
  
  # Go through all layer columns and replace all species codes with ftype
  for (i in 1:nrow(ftype)){
    for (layer in 1:n_layers){
      lpi[lpi[,layers[layer]] == ftype[i,1],layers[layer]] <- ftype[i,2]
    }
    # Replace codes in soilsurface
    if (ftype[i,1] %in% unique(lpi$soilsurface)){
      lpi[lpi$soilsurface == ftype[i,1], ]$soilsurface <- ftype[i,2]
    }
  }
  
  # Make lpi dataframe long
  lpi1 <- tidyr::pivot_longer(lpi[,c("plot","transect","point",layers, "soilsurface")], cols = names(lpi)[4:9], names_to = "layer")
  # Remove empty layers & layer column
  # Remove "S" and "N" as these will only be relevant for bareground calculation below
  lpi1 <- lpi1[!(lpi1$value %in% c("N","S","")),c(1:3,5)]
  
  # Now just get unique values
  lpi1 <- unique(lpi1)
  
  # Now, calculate number of times we see each species by plot
  plt_ct <- tapply(lpi1$value, lpi1$plot, function(x){table(x)})
  
  # Convert the list from above into a data frame
  nonzero_cover <- data.frame(plot = rep(names(plt_ct)[1], length(plt_ct[[1]])), 
                              species = plt_ct[1])
  names(nonzero_cover)[2:3] <- c("species","count")
  for (i in 2:length(plt_ct)){
    thisplot <- data.frame(plot = rep(names(plt_ct)[i], length(plt_ct[[i]])), 
                           species = plt_ct[i])
    names(thisplot)[2:3] <- c("species","count")
    nonzero_cover <- rbind(nonzero_cover,thisplot)
  }
  nonzero_cover$species <- as.character(nonzero_cover$species)
  
  # Special bareground calculation (does not consider litter as plant cover)
  # Should not consider rock pr other soilsurface as plant cover either
  NS <- lpi[lpi$toplayer == "N" & lpi$soilsurface %in% soilsurface,c(1,4:9)]
  bareground <- data.frame(table(NS$plot))
  bareground$species <- "bareground"
  bareground <- bareground[,c(1,3,2)]
  names(bareground) <- c("plot","species","count")
  
  # Add bareground to nonzero cover
  nonzero_cover <- rbind(nonzero_cover,bareground)
  
  # Calculate fractional cover
  nonzero_cover1 <- merge(nonzero_cover, pts_plt, by = "plot", all.x = T)
  nonzero_cover1$cover <- nonzero_cover1$count/nonzero_cover1$points
  
  # Reorder
  nonzero_cover1 <- nonzero_cover1[order(nonzero_cover1$plot, nonzero_cover1$species),]
  
  # Return columns of interest
  return(nonzero_cover1[,c(1,2,5)])
}