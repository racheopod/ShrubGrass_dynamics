#' Extract species-level cover from line-point intercept data.
#' 
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
#' @param separate_dead boolean. `TRUE` indicates that separate live and dead cover 
#'  values will be calculated for each species. `FALSE` indicates that a single
#'  cover value will be calculated for each species.
#' @param nonplants vector. A character vector of all non-plant non-soil surface
#'  codes present in data (default = c(`VL`,`HL`,`WL`)).
#' @param soilsurface vector. A character vector of all non-plant soil surface 
#'  codes present in the data (default = c(`BR`,`BY`,`CB`,`EL`,`GR`,`LC`,`M`,`S`,`ST`).
#' @param pres.abs_abund boolean. `TRUE` indicates a output will be a dataframe
#'  with rows corresponding to plots and columns corresponding to all species in 
#'  the dataset so that species not detected at a given plot in lpi will be 
#'  assigned a `0` (while fractional cover will be given for species present in 
#'  lpi). `FALSE` indicates that output will be a dataframe with three
#'  columns: `plot`,`species` (i.e. cover type), and `cover` (fractional). default
#'  is `FALSE`.
#' @return data frame of cover values for each plot. See `pre.abs_abund` for 
#'  details
#'  
#' @examples
#' # File path to final, cleaned data
#' cleandir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2022/02_Data_final'
#' 
#' # Read in lpi
#' lpi <- read.csv(file.path(cleandir,"lpi_2022.csv"))
#' 
#' #' Calculate cover (here we are just getting cover by species, not separating out live/dead)
#' cover <- lpicover(lpi, n_layers = 5, separate_dead = FALSE,
#'                   nonplants = c("VL","HL","WL"),
#'                   soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"),
#'                   pres.abs_abund = FALSE)
#' 
#' # Make a histogram of bareground, for example
#' par(mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
#' hist(cover[cover$species == "bareground",]$cover, main = "Histogram of bareground",
#'      xlab = "Bareground (fractional cover)", breaks = 10, xlim = c(0,1))
#' box()
#' abline(h = 0)


lpicover <- function(lpi, n_layers = 5, separate_dead = FALSE,
                     nonplants = c("VL","HL","WL"),
                     soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"),
                     pres.abs_abund = FALSE){
  # Check for soil surface plants that aren't found in layers above
  basalcheck <- lpi[!(lpi$soilsurface %in% soilsurface),]
  layers <- c("toplayer",paste0("layer",1:(n_layers-1)))
  if (sum(apply(basalcheck,1,function(x){!(x["soilsurface"] %in% x[layers])})) > 0){
    warning("One or more plant codes in soilsurface layer may not be present in upper layers.")
  }
  
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
  
  # Modify lpi data frame to be "long" instead of "wide"
  # First, in the case where 'separate_dead' == TRUE
  if (separate_dead){
    lpi_long <- as.matrix(lpi[,c("plot","toplayer", "toplayer_dead")])
    for (i in 1:(n_layers-1)){
      thislayer <- as.matrix(lpi[,c("plot",paste0("layer",i),paste0("layer",i,"_dead"))])
      lpi_long <- rbind(lpi_long, thislayer)
    }
    ssdf <- as.matrix(lpi[lpi$soilsurface %in% soilsurface,c("plot","soilsurface"),])
    ssdf <- cbind(ssdf,0)
    lpi_long <- rbind(lpi_long, ssdf)
    lpi_long <- as.data.frame(lpi_long)
    names(lpi_long)[2:3] <- c("layer","layer_dead")
    lpi_long$layer_dead <- as.numeric(lpi_long$layer_dead)
  } else if (!separate_dead){
    lpi_long <- as.matrix(lpi[,c("plot","toplayer")])
    for (i in 1:(n_layers-1)){
      thislayer <- as.matrix(lpi[,c("plot",paste0("layer",i))])
      lpi_long <- rbind(lpi_long, thislayer)
    }
    ssdf <- as.matrix(lpi[lpi$soilsurface %in% soilsurface,c("plot","soilsurface"),])
    lpi_long <- rbind(lpi_long, ssdf)
    lpi_long <- as.data.frame(lpi_long)
    names(lpi_long)[2] <- c("layer")
  }
  rownames(lpi_long) <- 1:nrow(lpi_long)
  
  # Remove empty observations (layers with no plant observed)
  lpi_long <- lpi_long[lpi_long$layer != "",]
  
  # Remove "S" and "N" as these will only be relevant for bareground calculation below
  lpi_long <- lpi_long[!(lpi_long$layer %in% c("N","S")),]
  
  # Calculate cover by species
  if (separate_dead){ # First when we want separate calculations for dead/live cover
    dead_status <- gsub("0","live",lpi_long$layer_dead)
    dead_status <- gsub("1","dead",dead_status)
    lpi_long$layer <- paste0(lpi_long$layer,"_",dead_status)
    lpi_long <- lpi_long[,-3]
    for (i in 1:length(nonplants)){
      lpi_long$layer <- gsub(paste0(nonplants[i],"_.*"),nonplants[i],lpi_long$layer)
    }
  }
  
  # Now, calculate number of times we see each species by plot
  plt_ct <- tapply(lpi_long$layer, lpi_long$plot, function(x){table(x)})
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

  # Special bareground calculation
  NS <- lpi[lpi$toplayer == "N" & lpi$soilsurface == "S" & 
              apply(lpi[,layers[-1]],1,function(x){unique(x)[1]  == ""}),c(1,4:9)]
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
  
  # Now build out a n_plots x n_species dataframe and populate with 0's where needed
  if (pres.abs_abund){
    species_df <- data.frame(plot = pts_plt$plot[1],species = sort(unique(nonzero_cover1$species)))
    for (i in 2:nrow(pts_plt)){
      thisdf <- data.frame(plot = pts_plt$plot[i],species = sort(unique(nonzero_cover1$species)))
      species_df <- rbind (species_df, thisdf)
    }
    pres.abs <- merge(species_df, nonzero_cover1[,c(1,2,5)], by = c("plot","species"), all.x = TRUE)
    pres.abs[is.na(pres.abs$cover),]$cover <- 0
    pres.abs1 <- reshape(pres.abs, timevar = "species", idvar = "plot", direction = "wide",
                         new.row.names = c(1:nrow(pts_plt)))
    names(pres.abs1) <- gsub("cover.","",names(pres.abs1))
    return(pres.abs1)
  } else{ # Or just return the simple cover x plot data frame for present species
    return(nonzero_cover1[,c(1,2,5)])
  }
}

#################################################
#' Make lpi data frame long for calculating cover (helper for `lpi_microsite` function)
#'
#' @param df lpi data frame with a `plot` column and columns for each layer
#' @param n_layers Number of layers.
#'
#' @return a long lpi data frame with all non-empty layers paired with plot name 
#'  and bound by row (`rbind()`). Output has column names `plot` and `layer`, with
#'  all records from all layers of each ftype or species in `layer` column.
#'

make_long <- function(df, n_layers){
  # Modify microsite data frames to be "long" instead of "wide"
  df_long <- as.matrix(df[,c("plot","toplayer")])
  for (i in 1:(n_layers-1)){
    thislayer <- as.matrix(df[,c("plot",paste0("layer",i))])
    df_long <- rbind(df_long, thislayer)
  }
  df_long <- as.data.frame(df_long)
  names(df_long)[2] <- c("layer")
  rownames(df_long) <- 1:nrow(df_long)

  # Remove empty observations (layers with no plant observed)
  return(df_long[df_long$layer != "",])
}


#################################################
#' Calculate cover from a long lpi data frame (helper for `lpi_microsite` function)
#'
#' @param df_long lpi data frame that has been converted to a "long" format using
#'   the `make_long` helper function above.
#' @param pts_plt Data frame with `plots`, `total_points` (total number of lpi hits
#'  for each plot), `canopy_points` (total number of points with shrub cover in 
#'  any layer), and `interspace_points` (total number of points without shrub cover).
#' 
#' @return A dataframe with columns for `plot`,`species`, `ms_count` (number of times
#'  this species was recorded in any layer), all columns from `pts_plt` (i.e., 
#'  `total_points`,`canopy_points`, and `interspace_points`), `total_cover` (`count`/
#'  `total_points`), and `ms_cover` (`count`/`canopy_points`).
#'

calc_cover <- function(df_long, pts_plt){
  # Now, calculate number of times we see each species by plot
  plt_ct <- tapply(df_long$layer, df_long$plot, function(x){table(x)})
  nonzero_cover <- data.frame(plot = rep(names(plt_ct)[1], length(plt_ct[[1]])), 
                              species = plt_ct[1])
  names(nonzero_cover)[2:3] <- c("species","ms_count")
  for (i in 2:length(plt_ct)){
    thisplot <- data.frame(plot = rep(names(plt_ct)[i], length(plt_ct[[i]])), 
                           species = plt_ct[i])
    names(thisplot)[2:3] <- c("species","ms_count")
    nonzero_cover <- rbind(nonzero_cover,thisplot)
  }
  nonzero_cover$species <- as.character(nonzero_cover$species)

  # Calculate fractional cover
  nonzero_cover1 <- merge(nonzero_cover, pts_plt, by = "plot", all.x = T)
  nonzero_cover1$total_cover <- nonzero_cover1$ms_count/nonzero_cover1$total_points
  nonzero_cover1$ms_cover <- nonzero_cover1$ms_count/nonzero_cover1$canopy_points
  return(nonzero_cover1)
}



#' Calculate herbaceous cover in shrub canopy and interspace microsites
#' 
#' @param lpi A dataframe of line-point intercept observations. This dataframe
#'  should include a "plot" column and columns named "toplayer","layer1","layer2",
#'  "layer3",... up to `n_layers`-1. If `separate_dead` is `TRUE`, there should 
#'  be additional columns corresponding to each layer column (named following the 
#'  pattern "toplayer_dead","layer1_dead","layer2_dead", etc.) that designate `1` 
#'  if the plant part in that layer was recorded as dead, and `0` otherwise.
#' @param n_layers Number of layers.
#' @param shrub_code The code in `ftype` that corresponds to a shrub. Default is
#'  `SH`.
#' @param ftypes Dataframe with a `species` column corresponding to all
#'  plant species codes in the `lpi` dataframe and a `ftype` column containing  
#'  the functional type of each species.#' 
#' @param soilsurface vector. A character vector of all non-plant soil surface 
#'  codes present in the data (default = c(`BR`,`BY`,`CB`,`EL`,`GR`,`LC`,`M`,`S`,
#'  `ST`).
#'  
#' @return A dataframe with columns for `plot`,`species`, `count` (number of times
#'  this species was recorded in any layer), `total_points` (total number of lpi
#'  points for each plot),`canopy_points` (total number of lpi points with shrub cover), 
#'  `interspace_points` (total number of lpi points without shrub cover), 
#'  `canopy_cover_of_total` (`canopy_count`/`total_points`), `canopy_cover` 
#'  (`canopy_count`/`canopy_points`), `interspace_cover_of_total` (`interspace_count`/
#'  `total_points`), `interspace_cover` (`interspace_count`/`interspace_points`),
#'  and `rii` (relative interaction index, calculated as (`canopy_cover` - 
#'  `interspace_cover`)/(`canopy_cover` + `interspace_cover`)).
#'  
#' @details Note that this returns a lot of information and in this case the 
#'  relative interaction index (`rii`) is calculated as the normalized ratio
#'  of the fractional covers of each plant type in canopy and interspace microsites.
#'  Negative values imply a competitive interaction (cover is higher in interspaces
#'  than in canopy microsites) while positive values imply a facultative interaction
#'  (cover is higher in canopy than in interspace microsites). This uses *relative* 
#'  cover within each microsite as opposed to absolute cover across the whole plot
#'  (which could be calculated manually using the `canopy_cover_of_total` and 
#'  `interspace_cover_of_total` variables if desired.)
#'  
#' @examples
#' # Set directory for cleaned data
#' cleandir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/2022/02_Data_final'
#' 
#' # Read in lpi
#' lpi <- read.csv(file.path(cleandir,"lpi_2022.csv"))
#'
#' # Create the ftype data frame
#' codes <- read.csv(file.path(cleandir,"codes_2022.csv"))
#' codes1 <- codes[,c("code","GrowthForm")]
#' names(codes1) <- c("species","ftype")
#' # Remove "S" (soil), "R" (rock), "N" (none), and "L" (litter)  rows
#' codes2 <- codes1[!(codes1$ftype %in% c("S","R","N","L")),]
#' ftype = codes2
#' 
#' # Calculate 
#' riis <- lpi_microsite(lpi, ftypes, n_layers = 5,
#'                       soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST"))

lpi_microsite <- function(lpi, ftypes, n_layers = 5, shrub_code = "SH",
                          soilsurface = c("BR","BY","CB","EL","GR","LC","M","S","ST")){
  # Check for soil surface plants that aren't found in layers above
  basalcheck <- lpi[!(lpi$soilsurface %in% soilsurface),]
  layers <- c("toplayer",paste0("layer",1:(n_layers-1)))
  if (sum(apply(basalcheck,1,function(x){!(x["soilsurface"] %in% x[layers])})) > 0){
    warning("One or more plant codes in soilsurface layer may not be present in upper layers.")
  }
  
  # Verify that there are no values == "" (empty) in top layer or soil surface
  if (sum(lpi$toplayer == "") > 0){
    stop(paste0("'toplayer' is missing for ",sum(lpi$toplayer == "")," records."))
  }
  if (sum(lpi$soilsurface == "") > 0){
    stop(paste0("'soilsurface' is missing for ",sum(lpi$soilsurface == "")," records."))
  }
  
  # Create points/plot named vector
  pts_plt1 <- data.frame(table(lpi$plot))
  names(pts_plt1) <- c("plot","total_points")
  
  # Remove records with "N" as top layer (no vegetation)
  foliar <- lpi[lpi$toplayer != "N",c("plot","transect","point",layers, "soilsurface")]
  
  # Go through all layer columns and replace all species codes with ftype
  for (i in 1:nrow(ftype)){
    for (layer in 1:n_layers){
      foliar[foliar[,layers[layer]] == ftype[i,1],layers[layer]] <- ftype[i,2]
      # Also remove any non-plant codes
      foliar[foliar[,layers[layer]] == "WL",layers[layer]] <- ""
      foliar[foliar[,layers[layer]] == "HL",layers[layer]] <- ""
      foliar[foliar[,layers[layer]] == "VL",layers[layer]] <- ""
    }
  }
  
  # Modify shrub_code for exact matches
  shrub_code1 <- paste0("\\b",shrub_code,"\\b")
  # Create a canopy column that indicates if there is a shrub anywhere in the hit
  foliar$canopy <- apply(foliar[,c(layers)], 1, function(x){as.numeric(sum(grepl(shrub_code1,x)) > 0)})
  
  # Separate canopy and interspace microsites
  canopy <- foliar[foliar$canopy == 1, ]
  interspace <- foliar[foliar$canopy == 0, ]
  
  # Calculate number of each kind of microsite per plot
  shrub_hits <- data.frame(table(canopy$plot))
  names(shrub_hits) <- c("plot","canopy_points")
  
  # Merge all points together
  pts_plt <- merge(pts_plt1, shrub_hits, by = "plot", all.x = TRUE)
  pts_plt$interspace_points <- pts_plt$total_points - pts_plt$canopy_points
  
  # Make canopy and interspace long, then calculate cover
  canopy1 <- calc_cover(make_long(canopy, n_layers), pts_plt)
  names(canopy1)[c(3,7:8)] <- c("canopy_count","canopy_cover_of_total","canopy_cover")
  interspace1 <- calc_cover(make_long(interspace, n_layers), pts_plt)
  names(interspace1)[c(3,7:8)] <- c("interspace_count","interspace_cover_of_total","interspace_cover")
  
  # Merge canopy and interspace
  cover <- merge(canopy1, interspace1[,c(1,2,3,7:8)], by = c("plot","species"), all = T)
  # Convert NA's in total_points to appropriate number
  for (i in 1:nrow(cover)){
    if (is.na(cover$total_points[i])){
      cover$total_points[i] <- pts_plt[pts_plt$plot == cover$plot[i],]$total_points
    }
  }
  # Convert NA's to 0 
  for (i in c(3,5:11)){
    cover[is.na(cover[,i]),i] <- 0
  }
  # Calculate RII
  cover$rii <- (cover$canopy_cover - cover$interspace_cover)/(cover$canopy_cover + cover$interspace_cover)
  return(cover)
}
