#' Calculate species or functional type average cover, height, and number for 
#' microsites (canopy and interspace) and canopy conditions (open and covered).
#' 
#' 
#' 
#' @param pcomm A dataframe of plant community observations. This dataframe
#'  should include columns for `plot`, `transect`, and `point`, followed by  
#'  columns ...
#' @param ftype
#' @param can_edge
#' @return data frame of average cover, height, and number for each species or 
#'  functional type category separated by microsite and canopy condition at each
#'  plot.
#' @details Average plot values will be computed across all observations recorded
#'  as each microsite or canopy condition, including points where there is a missing
#'  observation for one of the microsites (e.g., canopy microsite recorded but 
#'  interspace is missing for GP5-005 T3 25 m), or any point where there were no
#'  gaps > 60 cm (no interspace recorded). The number of quadrats used for each 
#'  calculated value is included in the output data frame for reference. Importantly,
#'  because the area surveyed in each canopy condition (0, 0.5, or 1 as proportion
#'  canopy cover over each observed plant) are most likely unequal, so we should
#'  only compare average cover/individual and height/individual, not *number* of 
#'  individual between canopy conditions.
#'  
#' @examples
#' # File path to final, cleaned data
#' datadir <- 'C:/Users/rache/Dropbox/Doctoral_projects/Projects/Shrub-grass_relationships/03_Data_Entry/Final_data'
#' 
#' # Read in plant community data
#' pcomm<- read.csv(file.path(datadir,"community.csv"))
#' 
#' # Create the ftype data frame
#' codes <- read.csv(file.path(datadir,"codes.csv"))
#' codes1 <- codes[,c("code","GrowthForm")]
#' names(codes1) <- c("species","ftype")
#' # Remove "S" (soil), "R" (rock), "N" (none), and "L" (litter)  rows
#' codes2 <- codes1[!(codes1$ftype %in% c("S","R","N","L")),]
#' ftype = codes2
#' 
#' # Calculate both microsite and canopy condition values
#' ms_and_ca_calcs <- pcomm_microsite(pcomm, ftype, can_edge = FALSE)
#' 
#' # Separate out each one
#' microsites <- ms_and_ca_calcs[[1]]
#' canopy_condition <- ms_and_ca_calcs[[2]]
#' 
#' 

pcomm_microsite <- function(pcomm, ftype, can_edge = FALSE){
  # Merge codes with community
  com <- merge(pcomm, ftype, by = "species", all.x = T)
  
  # Get unique ftypes
  ftypes <- sort(unique(com$ftype))
  ftypes <- ftypes[-which(ftypes %in% c("NO GAPS","SH"))]

  # Remove observations with no gaps (leave corresponding canopy cover microsite)
  com1 <- com[-which(!is.na(com$ftype) & com$ftype == "NO GAPS"), ]

  # Change from 0,1 to "Open", "Covered" if can_edge = TRUE
  if (can_edge){
    com1$canopy.condition <- NA
    com1[com1$canopy == 0,]$canopy.condition <- "open"
    com1[com1$canopy == 1,]$canopy.condition <- "covered"
    com1[com1$canopy == 0.5,]$canopy.condition <- "edge"
  } else {
    com1$canopy.condition <- NA
    com1[com1$canopy == 0,]$canopy.condition <- "open"
    com1[com1$canopy == 1 | com1$canopy == 0.5,]$canopy.condition <- "covered"
  }

  # Set number of plants to 0 for a "NO PLANTS" observation
  com1[com1$species == "NO PLANTS",]$number <- 0
  
  # Add cover column to com1
  com1$size <- com1$length * com1$width

  # Set size to 0.25 cm2 for rhizomatous grass
  # (if number = 0 it is a "NO PLANTS" observation)
  for (i in 1:nrow(com1)){
    if (com1$number[i] != 0 && is.na(com1$size[i])){
      com1$size[i] <- 0.25 
    } 
    # Set size to 0.25 cm2 for rhizomatous annual plants
    # And make sure length & width is NA to avoid confusion
    if (!is.na(com1$ftype[i]) && (com1$ftype[i] == "AF" | com1$ftype[i] == "AG")){
      com1$size[i] <- 0.25
      com1$length[i] <- NA
      com1$width[i] <- NA
    }
    # And make size 0 when there were no plants observed in that microsite
    if (com1$species[i] == "NO PLANTS"){
      com1$size[i] <- 0
    } 
  }
  
  # Loop through ftypes and calculate avg cover, heigh, number/ microsite and
  # avg cover, height / canopy condition
  results_list <-list()
  results_list[[1]] <- pft_summary(thispft = com1[!is.na(com1$ftype) & com1$ftype == ftypes[1],], 
                         pcomm, thisft = ftypes[1])
  for (ft in 2:length(ftypes)){
    results_list[[ft]] <- pft_summary(thispft = com1[!is.na(com1$ftype) & com1$ftype == ftypes[ft],], 
                                      pcomm, thisft = ftypes[ft])
  }
  
  # Unlist and make two dfs to return
  microsites <- results_list[[1]][[1]]
  canopy_condition <- results_list[[1]][[2]]
  for (i in 2:length(results_list)){
    microsites <- rbind(microsites, results_list[[i]][[1]])
    canopy_condition <- rbind(canopy_condition, results_list[[i]][[2]])
  }
  return(list(microsites, canopy_condition))
}

 
#' Helper function for `pcomm_microsite`, calculates average cover, height, and 
#' number as appropriate for microsites and canopy conditions
#' 
#' 
#' 
#' @param thispft A dataframe that contains only `pcomm` observations for a single
#'  species or functional type.
#' @param pcomm The full `pcomm` data frame.
#' @param thisft The functional type (or species) to use for labeling the calculated
#'  values in the output data frames.
#' @return A list of two data frames containing information for microsites and canopy
#'  conditions for `thisft`. The first is the `microsite` data frame and the second
#'  is the `canopy_condition` data frame.
#' 

pft_summary <- function(thispft, pcomm, thisft = ftypes[ft]){
  # Pull in pipe from magrittr
  `%>%` <- magrittr::`%>%`
  # Make plots data frame for all microsites that were sampled:
  plots_ms <- unique(pcomm[,c("plot","transect","point","microsite")])
  # Make plots data frame for all plots
  plots_ca <- data.frame(unique(pcomm$plot))
  
  # First summarize by microsite
  ms <- thispft %>%
        dplyr::group_by(plot, transect, point, microsite) %>%
        dplyr::summarize(number = sum(number),
                  cover = sum(size)/(20*30),
                  size = sum(size), # Add together all sizes to average across number in next step
                  height = sum(height)) # Same as size
  ms$code <- NA
  # Add in the most common species for each microsite
  for (r in 1:nrow(ms)){
     ms$code[r] <- names(table(thispft[thispft$plot == ms$plot[r] & thispft$transect == ms$transect[r] &
                            thispft$point == ms$point[r] & thispft$microsite == ms$microsite[r],]$species)
            [which.max(table(thispft[thispft$plot == ms$plot[r] & thispft$transect == ms$transect[r] &
                      thispft$point == ms$point[r] & thispft$microsite == ms$microsite[r],]$species))])
  }
  
  # Merge with "plots" to account for cover and number where no plants were recorded
  ms <- merge(plots_ms, ms, by = c("plot","transect","point","microsite"), all = TRUE)
  # Fill zeros where needed (don't set height to 0 bc we don't want to average w/ 0)
  if (nrow(ms[is.na(ms$number),]) != 0){
    ms[is.na(ms$number),]$number <- 0
  }
  if (nrow(ms[is.na(ms$cover),]) != 0){
    ms[is.na(ms$cover),]$cover <- 0
  }

  # Now summarize just by plot
  ms1 <- ms %>%
    dplyr::group_by(plot, microsite) %>%
    dplyr::summarize(n_quadrats = dplyr::n(),
              avg_number = mean(number),
              avg_cover = mean(cover),
              n_ind = sum(number),
              size = sum(size, na.rm = TRUE)/sum(number),
              height = sum(height, na.rm = TRUE)/sum(number))
  ms1$code <- NA
  # Add in the most common species for each microsite
  for (r in 1:nrow(ms1)){
    if (dim(table(thispft[thispft$plot == ms1$plot[r] & thispft$microsite == ms1$microsite[r],]$species))>0){
      ms1$code[r] <- names(table(thispft[thispft$plot == ms1$plot[r] & thispft$microsite == ms1$microsite[r],]$species)
                          [which.max(table(thispft[thispft$plot == ms1$plot[r] & thispft$microsite == ms1$microsite[r],]$species))])
    } else{
      ms1$code[r] <- NA
    }
  }
  ms1 <- as.data.frame(ms1)
  
  #########################################
  # Then summarize by canopy condition
  ca <- thispft %>%
    dplyr::group_by(plot, canopy.condition) %>%
    dplyr::summarize(n_ind = sum(number),
              avg_size = mean(size),
              height = mean(height))
  ca$code <- NA
  # Add in the most common species for each microsite
  for (r in 1:nrow(ca)){
    ca$code[r] <- names(table(thispft[thispft$plot == ca$plot[r] & thispft$canopy.condition == ca$canopy.condition[r],]$species)
                        [which.max(table(thispft[thispft$plot == ca$plot[r] & thispft$canopy.condition == ca$canopy.condition[r],]$species))])
  }
  ca <- as.data.frame(ca)
  
  #########################################
  # Reshape dfs
  ms2 <- reshape(ms1, direction="wide", idvar = "plot", timevar = "microsite")
  ms2$ftype <- thisft
  names(ms2)[2:15] <- c("canopy_n_quadrats","canopy_avg_number","canopy_avg_cover",
                        "canopy_n_ind","canopy_size","canopy_height","canopy_code",
                       "interspace_n_quadrats","interspace_avg_number","interspace_avg_cover",
                       "interspace_n_ind","interspace_size","interspace_height",
                       "interspace_code")
  ms3 <- ms2[,c("plot","ftype","canopy_n_quadrats","interspace_n_quadrats",
                "canopy_avg_number","interspace_avg_number",
                "canopy_avg_cover","interspace_avg_cover",
                "canopy_n_ind","interspace_n_ind","canopy_size","interspace_size",
                "canopy_height","interspace_height","canopy_code","interspace_code")]
  
  ca1 <- reshape(ca, direction="wide", idvar = "plot", timevar = "canopy.condition")
  ca1$ftype <- thisft
  names(ca1)[2:9] <- c("open_n_ind","open_size","open_height","open_code",
                       "canopy_n_ind","canopy_size","canopy_height","canopy_code")
  ca2 <- ca1[,c("plot","ftype","canopy_n_ind","open_n_ind","canopy_size","open_size",
                "canopy_height","open_height","open_code","canopy_code")]
  
  return(list(ms3,ca2))
}

#' Calculate relative interaction index for any pair of values from canopy and
#' interspace microsites
#' 
#' @param canopy Vector of values corresponding to plants in canopy microsites
#'  (or in the 'canopy' condition, i.e., shrub canopy cover over plant).
#' @param interspace Vector of values corresponding to plants in interspace microsites
#'  (or in the 'open' condition, i.e., no shrub canopy over plant)
#' @return Vector of the relative interaction index (rii) calculated from each 
#'  pair of `canopy` and `interspace` values. RII = (canopy - interspace)/
#'  (canopy + interspace), and ranges from -1 (strong competition) to 1 (strong
#'  facilitation). *Note of caution:* When the values in `canopy` and `interspace` 
#'  are both 0, the RII = 1 and will suggest strong facilitation, but the nature 
#'  of the interaction of shrubs with 0 plants cannot be facultative--so watch out.
#' 

calculate_rii <- function(canopy, interspace){
    return((canopy-interspace)/(canopy+interspace))
}