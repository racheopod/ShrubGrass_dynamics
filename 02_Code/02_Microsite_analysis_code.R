#####################################
# Purpose: Analyze abundance/size of different plant types in microsites from Summer 2021 data
#
#
# Rachel R. Renne
# March 9, 2022
#
##################################

library(dplyr)

# Read in community data
community <- read.csv("C:/Users/rache/Dropbox/Doctoral Projects/Data/Final Data 2022/community.csv")
# Remove GP5-005 I obs from T3, point 25 (due to missing Canopy microsite)
community <- community[-which(community$plot == "GP5-005" & community$transect == 3 & community$point == 25),]

# Read in plot data (with MAP)
plots <- read.csv("C:/Users/rache/Dropbox/Doctoral Projects/Data/Final Data 2022/plots.csv")

# Read in plant list
plants <- read.csv("C:/Users/rache/Dropbox/Doctoral Projects/Data/Final Data 2022/codes_edited.csv")
# Get just codes
plants1 <- plants[,c(1,5:6)]

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

# Convert all 0.5 canopy conditions to 1
com1[com1$canopy == 0.5,]$canopy <- 1

# Change from 0,1 to "Open", "Covered"
com1[com1$canopy == 0,]$canopy <- "O"
com1[com1$canopy == 1,]$canopy <- "C"


# Troubleshoot unique plot-transect-point combinations
unq <- unique(com1[,c(2,3,6)])
unq <- unq[order(unq$plot, unq$transect, unq$point),]
table(unq$plot, unq$transect)
# OK, we are missing point 25 on transect 2 for WYBas3-014 (and removed GP5-005 T3 P25)

# Create automated PFT analysis function
ft_analysis <- function(pfts, unq, pftname){
  
  # Make a results data frame to hold bootstrap samples
  results <- data.frame(plot=plots$plot, ms_diff_numplants = NA, ca_diff_numplants = NA,
                        ms_diff_avgheight = NA, ca_diff_avgheight = NA,
                        ms_diff_cover = NA, ca_diff_cover = NA)
  
  # Make another results df for storing results from lm's
  results1 <- data.frame(pft = rep(pftname, 6), 
                         variable = c("numplants","numplants",
                                      "avgheight","avgheight","cover","cover"), 
                         factor = rep(c("ms","ca"),3), p_map = NA, 
                         p_factor = NA, p_interaction = NA)
  
  # Pull out this/these PFTs
  thispft <- com1[com1$GrowthForm %in% pfts | com1$GrowthForm == "NO PLANTS",]

  # First summarize by microsite
  ms <- thispft %>%
    group_by(plot, transect, point, microsite) %>%
    summarize(num_plants = sum(number, na.rm = T),
              avg_height = mean(height, na.rm = T),
              cover = sum(length*width),
              code = names(table(Code)[which.max(table(Code))]))
  # Add back in microsites that were lost
  all_plots <- cbind(rbind(unq, unq), c(rep("C",nrow(unq)),rep("I",nrow(unq))))
  names(all_plots)[4] <- "microsite"
  ms <- merge(ms, all_plots, by = c("plot","transect","point","microsite"), all = TRUE)

  # Convert number of plants to 0 where no observations were recorded for PG
  if (nrow(ms[is.na(ms$num_plants),]) != 0){
  ms[is.na(ms$num_plants),]$num_plants <- 0
  }

  # Make small cover amounts (0.25 cm2) for rhizomatous grasses
  # And make cover 0 when there were no plants observed in that microsite
  for (i in 1:nrow(ms)){
    if (ms$num_plants[i] != 0 && is.na(ms$cover[i])){
      ms$cover[i] <- 0.25 
    } 
    if (is.na(ms$code[i]) | ms$code[i] == "NO PLANTS"){
      ms$cover[i] <- 0
    } 
  }

  # Then summarize by canopy condition
  ca <- thispft %>%
    group_by(plot, transect, point, canopy) %>%
    summarize(num_plants = sum(number, na.rm = T),
              avg_height = mean(height, na.rm = T),
              cover = sum(length*width),
              code = names(table(Code)[which.max(table(Code))]))
  all_plots <- cbind(rbind(unq, unq),c(rep("O",nrow(unq)),rep("C",nrow(unq))))
  names(all_plots)[4] <- "canopy"
  ca <- merge(ca, all_plots, by = c("plot","transect","point","canopy"), all = TRUE)

  # Convert number of plants to 0 where no observations were recorded for PG
  if (nrow(ca[is.na(ca$num_plants),]) != 0){
  ca[is.na(ca$num_plants),]$num_plants <- 0
  }
  
  # Make small cover amounts (0.25 cm2) for rhizomatous grasses
  # And make cover 0 when there were no plants observed in that microsite
  for (i in 1:nrow(ca)){
    if (ca$num_plants[i] != 0 && is.na(ca$cover[i])){
      ca$cover[i] <- 0.25 
    } 
    if (is.na(ca$code[i]) | ca$code[i] == "NO PLANTS"){
      ca$cover[i] <- 0
    } 
  }

  # Finally, merge with plots to get MAP
  ms <- merge(ms, plots[,c("plot","altitude","map","mat")], by = "plot", all.x = T)
  ca <- merge(ca, plots[,c("plot","altitude","map","mat")], by = "plot", all.x = T)


  ################################################
  # A. Number of plants
  ################################################
  
  print("Now working on number of plants.")
  
  #################################
  # Microsite
  #################################
  # T.test overall:
  ms_ttest <- t.test(ms$num_plants~ms$microsite)

  plots_numplants_ms <- NULL
  for (i in 1:17){
    thisplot <- ms[ms$plot == plots$plot[i],]
    diffx <- sum(thisplot[thisplot$microsite == "C", ]$num_plants) -
             sum(thisplot[thisplot$microsite == "I", ]$num_plants)
    bs <- NULL
    for (b in 1:10000){
      thisplot$microsite <- sample(thisplot$microsite)
      bs[b] <- sum(thisplot[thisplot$microsite == "C", ]$num_plants) -
               sum(thisplot[thisplot$microsite == "I", ]$num_plants)
    }
    plots_numplants_ms[i] <- mean(abs(bs) >= abs(diffx), na.rm = T)
  }

  results$ms_diff_numplants <- plots_numplants_ms

  #################################
  # Canopy cover
  #################################
  # T.test overall:
  ca_ttest <- t.test(ca$num_plants~ca$canopy)

  plots_numplants_ca <- NULL
  for (i in 1:17){
    thisplot <- ca[ca$plot == plots$plot[i],]
    diffx <- sum(thisplot[thisplot$canopy == "C", ]$num_plants) -
      sum(thisplot[thisplot$canopy == "O", ]$num_plants)
    bs <- NULL
    for (b in 1:10000){
      thisplot$canopy <- sample(thisplot$canopy)
      bs[b] <- sum(thisplot[thisplot$canopy == "C", ]$num_plants) -
      sum(thisplot[thisplot$canopy == "O", ]$num_plants)
    }
  plots_numplants_ca[i] <- mean(abs(bs) >= abs(diffx), na.rm = T)
  }

  results$ca_diff_numplants <- plots_numplants_ca

  png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_num_plants.png"), 
      width = 16, height = 16, units = "in", res = 300)
  par(mfrow = c(4,10))
  for (i in 1:17){
    ms1 <- ms[ms$plot == plots$plot[i],]
    ca1 <- ca[ca$plot == plots$plot[i],]
    par(mar = c(2,2.2,2,0), tcl = 0.3, mgp = c(1,0,0))
    boxplot(ms1$num_plants~ms1$microsite, cex.lab = 1.5,
            main = "Microsite",xlab = "",ylab = "Number of plants", cex.axis = 1,
            names = c("Canopy", "Interspace"))
    legend("topleft", legend = c(plots$plot[i], paste0("p = ",round(results[i,2],2))), bty = "n")
    par(mar = c(2,0,2,2.2))
    boxplot(ca1$num_plants~ca1$canopy, cex.lab = 1.5,
            main = "Canopy cover",xlab = "",yaxt = "n", cex.axis=1,
            names = c("Covered", "Open"))
    legend("topleft", legend = c("", paste0("p = ",round(results[i,3],2))), bty = "n")
    
  }
  dev.off()

  png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_num_plants_combined.png"),
      width = 4, height = 4, units = "in", res = 300)
  par(mfrow = c(1,2), mar = c(2,2.2,2,0), tcl = 0.3, mgp = c(1,0,0))
  boxplot(ms$num_plants~ms$microsite, cex.lab = 1.5,
          main = "Microsite",xlab = "",ylab = "Number of plants", cex.axis = 0.8,
          names = c("Canopy", "Interspace"),ylim = c(0,40))
  legend("topleft", legend = paste0("p = ",round(ms_ttest$p.value,2)), bty = "n")
  par(mar = c(2,0,2,2.2))
  boxplot(ca$num_plants~ca$canopy, cex.lab = 1.5, names = c("Covered","Open"),
          main = "Canopy cover",xlab = "",ylim = c(0,40),yaxt = "n", cex.axis=0.8)
  legend("topleft", legend = paste0("p = ",round(ca_ttest$p.value,2)), bty = "n")
  dev.off()

  ##################
  # Model with precipitation
  ##################

  ms_num <- summary(lm(ms$num_plants~ms$map*ms$microsite))
  ca_num <- summary(lm(ca$num_plants~ca$map*ca$canopy))
  
  #Save results
  results1[1,4:6] <- ms_num$`Pr(>F)`[1:3]
  results1[2,4:6] <- ca_num$`Pr(>F)`[1:3]

  ################################################
  # B. Height
  ################################################

  print("Now working on average height.")
  
  #################################
  # Microsite
  #################################

  # T.test
  ms_ttest <- t.test(ms$avg_height~ms$microsite)

  plots_height_ms <- NULL
  for (i in 1:17){
    thisplot <- ms[ms$plot == plots$plot[i],]
    diffx <- median(thisplot[thisplot$microsite == "C", ]$avg_height, na.rm = T) -
      median(thisplot[thisplot$microsite == "I", ]$avg_height, na.rm = T)
    bs <- NULL
    for (b in 1:10000){
      thisplot$microsite <- sample(thisplot$microsite)
      bs[b] <- median(thisplot[thisplot$microsite == "C", ]$avg_height, na.rm = T) -
        median(thisplot[thisplot$microsite == "I", ]$avg_height, na.rm = T)
    }
    plots_height_ms[i] <- mean(abs(bs) >= abs(diffx), na.rm = T)
  }

  results$ms_diff_avgheight <- plots_height_ms

  #################################
  # Canopy cover
  #################################
  # T.test overall:
  ca_ttest <- t.test(ca$avg_height~ca$canopy)

  plots_height_ca <- NULL
  for (i in 1:17){
    thisplot <- ca[ca$plot == plots$plot[i],]
    diffx <- median(thisplot[thisplot$canopy == "O", ]$avg_height, na.rm = T) -
      median(thisplot[thisplot$canopy == "C", ]$avg_height, na.rm = T)
    bs <- NULL
    for (b in 1:10000){
      thisplot$canopy <- sample(thisplot$canopy)
      bs[b] <- median(thisplot[thisplot$canopy == "O", ]$avg_height, na.rm = T) -
        median(thisplot[thisplot$canopy == "C", ]$avg_height, na.rm = T)
    }
    plots_height_ca[i] <- mean(abs(bs) >= abs(diffx), na.rm = T)
  }

  results$ca_diff_avgheight <- plots_height_ca

  png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_avgheight.png"), 
      width = 16, height = 16, units = "in", res = 300)
  par(mfrow = c(4,10))
  for (i in 1:17){
    ms1 <- ms[ms$plot == plots$plot[i],]
    ca1 <- ca[ca$plot == plots$plot[i],]
    par(mar = c(2,2.2,2,0), tcl = 0.3, mgp = c(1,0,0))
    if (length(table(ms1[!is.na(ms1$avg_height),]$microsite)) == 2){
    boxplot(ms1$avg_height~ms1$microsite, cex.lab = 1.5,
            main = "Microsite",xlab = "",ylab = "Height (cm)", cex.axis = 1,
            names = c("Canopy", "Interspace"),ylim = c(0,45))
    } else{
      boxplot(c(1:10) ~ c(rep(1,5),rep(2,5)), col = "white", cex.lab = 1.5,
              border = "white", fill = "white",
           main = "Microsite",xlab = "",ylab = "Height (cm)", cex.axis = 1,
           names = c("Canopy", "Interspace"),ylim = c(0,45))
      legend("center", legend = "NO PLANTS", bty = "n")
    }
    legend("topleft", legend = c(plots$plot[i], paste0("p = ",round(results[i,4],2))), bty = "n")
    par(mar = c(2,0,2,2.2))
    if (length(table(ca1[!is.na(ca1$avg_height),]$canopy)) == 2){
    boxplot(ca1$avg_height~ca1$canopy, cex.lab = 1.5, names = c("Covered","Open"),
          main = "Canopy cover",xlab = "",ylim = c(0,45),yaxt = "n", cex.axis=1)
    } else {
      boxplot(c(ca1$avg_height,0)~c(ca1$canopy,"O"), cex.lab = 1.5, names = c("Covered","Open"),
              main = "Canopy cover",xlab = "",ylim = c(0,45),yaxt = "n", cex.axis=1,
              border = c("black","white"))
      legend("bottomright", legend = "NO PLANTS", bty = "n")
    }
    legend("topleft", legend = c("", paste0("p = ",round(results[i,5],2))), bty = "n")
    
  }
  dev.off()
  
  png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_avgheight_combined.png"),
      width = 4, height = 4, units = "in", res = 300)
  par(mar = c(2,2.2,2,0), tcl = 0.3, mgp = c(1,0,0), mfrow = c(1,2))
  boxplot(ms$avg_height~ms$microsite, cex.lab = 1.5,
          main = "Microsite",xlab = "",ylab = "Height (cm)", cex.axis = 0.8,
          names = c("Canopy", "Interspace"),ylim = c(0,45))
  legend("topleft", legend = paste0("p = ",round(ms_ttest$p.value,2)), bty = "n")
  par(mar = c(2,0,2,2.2))
  boxplot(ca$avg_height~ca$canopy, cex.lab = 1.5,names=c("Covered","Open"),
          main = "Canopy cover",xlab = "",ylim = c(0,45),yaxt = "n", cex.axis=0.8)
  legend("topleft", legend = paste0("p = ",round(ca_ttest$p.value,2)), bty = "n")
  dev.off()
  
  ##################
  # Model with precipitation
  ##################
  
  ms_ht <- summary(lm(ms$avg_height~ms$map*ms$microsite))
  ca_ht <- summary(lm(ca$avg_height~ca$map*ca$canopy))
  
  #Save results
  results1[3,4:6] <- ms_ht$`Pr(>F)`[1:3]
  results1[4,4:6] <- ca_ht$`Pr(>F)`[1:3]
  


  ################################################
  # C. Cover
  ################################################
  
  print("Now working on cover.")

  #################################
  # Microsite
  #################################

  # T.test
  ms_ttest <- t.test(ms$cover~ms$microsite)

  plots_cover_ms <- NULL
  for (i in 1:17){
    thisplot <- ms[ms$plot == plots$plot[i],]
    diffx <- median(thisplot[thisplot$microsite == "C", ]$cover, na.rm = T) -
      median(thisplot[thisplot$microsite == "I", ]$cover, na.rm = T)
    bs <- NULL
    for (b in 1:10000){
      thisplot$microsite <- sample(thisplot$microsite)
      bs[b] <- median(thisplot[thisplot$microsite == "C", ]$cover, na.rm = T) -
        median(thisplot[thisplot$microsite == "I", ]$cover, na.rm = T)
    }
    plots_cover_ms[i] <- mean(abs(bs) >= abs(diffx), na.rm = T)
  }

  results$ms_diff_cover <- plots_cover_ms
  
  
  #################################
  # Canopy cover
  #################################
  # T.test overall:
  ca_ttest <- t.test(ca$cover~ca$canopy)

  plots_cover_ca <- NULL
  for (i in 1:17){
    thisplot <- ca[ca$plot == plots$plot[i],]
    diffx <- median(thisplot[thisplot$canopy == "O", ]$cover, na.rm = T) -
      median(thisplot[thisplot$canopy == "C", ]$cover, na.rm = T)
    bs <- NULL
    for (b in 1:10000){
      thisplot$canopy <- sample(thisplot$canopy)
      bs[b] <- median(thisplot[thisplot$canopy == "O", ]$cover, na.rm = T) -
        median(thisplot[thisplot$canopy == "C", ]$cover, na.rm = T)
    }
    plots_cover_ca[i] <- mean(abs(bs) >= abs(diffx), na.rm = T)
  }


  results$ca_diff_cover <- plots_cover_ca

  
  png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_cover.png"), 
      width = 16, height = 16, units = "in", res = 300)
  par(mfrow = c(4,10))
  for (i in 1:17){
    ms1 <- ms[ms$plot == plots$plot[i],]
    ca1 <- ca[ca$plot == plots$plot[i],]
    par(mar = c(2,2.2,2,0), tcl = 0.3, mgp = c(1,0,0))
    boxplot(ms1$cover~ms1$microsite, cex.lab = 1.5,
            main = "Microsite",xlab = "",ylab = "Cover (cm2)", cex.axis = 1,
            names = c("Canopy", "Interspace"))
    legend("topleft", legend = c(plots$plot[i], paste0("p = ",round(results[i,6],2))), bty = "n")
    par(mar = c(2,0,2,2.2))
    boxplot(ca1$cover~ca1$canopy, cex.lab = 1.5,names=c("Covered","Open"),
            main = "Canopy cover",xlab = "",yaxt = "n", cex.axis=1)
    legend("topleft", legend = c("", paste0("p = ",round(results[i,7],2))), bty = "n")
  }
  dev.off()

  png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_cover_combined.png"),
      width = 4, height = 4, units = "in", res = 300)
  par(mfrow = c(1,2))
  par(mar = c(2,2.2,2,0), tcl = 0.3, mgp = c(1,0,0))
  boxplot(ms$cover~ms$microsite, cex.lab = 1.5,
          main = "Microsite",xlab = "",ylab = "Cover (cm2)", cex.axis = 0.8,
          names = c("Canopy", "Interspace"), ylim = c(0,325))
  legend("topleft", legend = paste0("p = ",round(ms_ttest$p.value,2)), bty = "n")
  par(mar = c(2,0,2,2.2))
  boxplot(ca$cover~ca$canopy, cex.lab = 1.5,ylim=c(0,325),
          names = c("Covered","Open"),
          main = "Canopy cover",xlab = "",yaxt = "n", cex.axis=0.8)
  legend("topleft", legend = paste0("p = ",round(ca_ttest$p.value,2)), bty = "n")
  dev.off()
  
  ##################
  # Model with precipitation
  ##################
  
  ms_cov <- summary(lm(ms$cover~ms$map*ms$microsite))
  ca_cov <- summary(lm(ca$cover~ca$map*ca$canopy))
  
  #Save results
  results1[5,4:6] <- ms_cov$`Pr(>F)`[1:3]
  results1[6,4:6] <- ca_cov$`Pr(>F)`[1:3]
 
  # Save results to file
  write.csv(results, paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_bootstrap.csv"),
            row.names = FALSE)
  
  return(results1)
}
  
# Run on different pfts

# All plants:
ft_analysis(pfts = c("PG","PF","AF","SSH","AG","SE","SUC"),unq, pftname = "allplants")
  
# Pherbs
ft_analysis(pfts = c("PG","PF","SE"), unq, pftname = "Pherbs")
  
# Perennial grasses
ft_analysis(pfts=c("PG"), unq, pftname = "PG")

# Perennial forbs
ft_analysis(pfts=c("PF"), unq, pftname = "PF")
  
  
#################################
# Make some models--lm
#################################


pickmodel <- function(model_list, x){
  if (anova(model_list[[4]])$`Pr(>F)`[3] < 0.05){
    y = model_list[[4]]$coefficients[1]+model_list[[4]]$coefficients[2]*x
    lines(y~x, col = rgb(0,0.5,0,1), lwd = 2)
    y = model_list[[4]]$coefficients[1]+model_list[[4]]$coefficients[3]+
      model_list[[4]]$coefficients[2]*x+model_list[[4]]$coefficients[4]*x
    lines(y~x, col = rgb(0.6,0.3,0.2,1), lwd = 2)
  }else if (anova(model_list[[2]])$`Pr(>F)`[1] < 0.05 && anova(model_list[[1]])$`Pr(>F)`[1] < 0.05){
    y = model_list[[3]]$coefficients[1]+model_list[[3]]$coefficients[2]*x
    lines(y~x, col = rgb(0,0.5,0,1), lwd = 2)
    y = model_list[[3]]$coefficients[1]+model_list[[3]]$coefficients[3]+
      model_list[[3]]$coefficients[2]*x
    lines(y~x, col = rgb(0.6,0.3,0.2,1), lwd = 2)
  } else if (anova(model_list[[2]])$`Pr(>F)`[1] < 0.05 && anova(model_list[[1]])$`Pr(>F)`[1] > 0.05){
    lines(rep(model_list[[2]]$coefficients[1],length(x))~x, col = rgb(0,0.5,0,1), lwd = 2)
    lines(rep(model_list[[2]]$coefficients[1] + model_list[[2]]$coefficients[2], length(x))~x, col = rgb(0.6,0.3,0.2,1), lwd = 2)
  } else if (anova(model_list[[2]])$`Pr(>F)`[1] > 0.05 && anova(model_list[[1]])$`Pr(>F)`[1] < 0.05){
    lines((model_list[[1]]$coefficients[2]*x + model_list[[1]]$coefficients[1])~x, lwd =2)
  }
}


# function to show linear models
plot_results <- function(pfts, x, pftname){

# Pull out this/these PFTs
thispft <- com1[com1$Code %in% pfts | com1$GrowthForm == "NO PLANTS",]

# First summarize by microsite
ms <- thispft %>%
  group_by(plot, transect, point, microsite) %>%
  summarize(num_plants = sum(number, na.rm = T),
            avg_height = mean(height, na.rm = T),
            cover = sum(length*width),
            code = names(table(Code)[which.max(table(Code))]))
# Add back in microsites that were lost
all_plots <- cbind(rbind(unq, unq), c(rep("C",nrow(unq)),rep("I",nrow(unq))))
names(all_plots)[4] <- "microsite"
ms <- merge(ms, all_plots, by = c("plot","transect","point","microsite"), all = TRUE)

# Convert number of plants to 0 where no observations were recorded for PG
if (nrow(ms[is.na(ms$num_plants),]) != 0){
  ms[is.na(ms$num_plants),]$num_plants <- 0
}

# Make small cover amounts (0.25 cm2) for rhizomatous grasses
# And make cover 0 when there were no plants observed in that microsite
for (i in 1:nrow(ms)){
  if (ms$num_plants[i] != 0 && is.na(ms$cover[i])){
    ms$cover[i] <- 0.25 
  } 
  if (is.na(ms$code[i]) | ms$code[i] == "NO PLANTS"){
    ms$cover[i] <- 0
  } 
}

# Then summarize by canopy condition
ca <- thispft %>%
  group_by(plot, transect, point, canopy) %>%
  summarize(num_plants = sum(number, na.rm = T),
            avg_height = mean(height, na.rm = T),
            cover = sum(length*width),
            code = names(table(Code)[which.max(table(Code))]))
all_plots <- cbind(rbind(unq, unq),c(rep("O",nrow(unq)),rep("C",nrow(unq))))
names(all_plots)[4] <- "canopy"
ca <- merge(ca, all_plots, by = c("plot","transect","point","canopy"), all = TRUE)

# Convert number of plants to 0 where no observations were recorded for PG
if (nrow(ca[is.na(ca$num_plants),]) != 0){
  ca[is.na(ca$num_plants),]$num_plants <- 0
}

# Make small cover amounts (0.25 cm2) for rhizomatous grasses
# And make cover 0 when there were no plants observed in that microsite
for (i in 1:nrow(ca)){
  if (ca$num_plants[i] != 0 && is.na(ca$cover[i])){
    ca$cover[i] <- 0.25 
  } 
  if (is.na(ca$code[i]) | ca$code[i] == "NO PLANTS"){
    ca$cover[i] <- 0
  } 
}

# Finally, merge with plots to get MAP
ms <- merge(ms, plots[,c("plot","altitude","map","mat")], by = "plot", all.x = T)
ca <- merge(ca, plots[,c("plot","altitude","map","mat")], by = "plot", all.x = T)


m_cov_map <- (lm(sqrt(ms$cover)~ms$map))
m_cov_ft <-  (lm(sqrt(ms$cover)~ms$microsite))
m_cov_add <- (lm(sqrt(ms$cover)~ms$map+ms$microsite)) 
m_cov_int <- (lm(sqrt(ms$cover)~ms$map*ms$microsite)) 
m_ht_map <- (lm(sqrt(ms$avg_height)~ms$map))
m_ht_ft <- (lm(sqrt(ms$avg_height)~ms$microsite))
m_ht_add <- (lm(sqrt(ms$avg_height)~ms$map+ms$microsite))
m_ht_int <- (lm(sqrt(ms$avg_height)~ms$map*ms$microsite))
m_num_map <- (lm(sqrt(ms$num_plants)~ms$map))
m_num_ft <- (lm(sqrt(ms$num_plants)~ms$microsite))
m_num_add <- (lm(sqrt(ms$num_plants)~ms$map+ms$microsite))
m_num_int <- (lm(sqrt(ms$num_plants)~ms$map*ms$microsite))

c_cov_map <- (lm(sqrt(ca$cover)~ca$map))
c_cov_ft <-  (lm(sqrt(ca$cover)~ca$canopy))
c_cov_add <- (lm(sqrt(ca$cover)~ca$map+ca$canopy)) 
c_cov_int <- (lm(sqrt(ca$cover)~ca$map*ca$canopy)) 
c_ht_map <- (lm(sqrt(ca$avg_height)~ca$map))
c_ht_ft <- (lm(sqrt(ca$avg_height)~ca$canopy))
c_ht_add <- (lm(sqrt(ca$avg_height)~ca$map+ca$canopy))
c_ht_int <- (lm(sqrt(ca$avg_height)~ca$map*ca$canopy))
c_num_map <- (lm(sqrt(ca$num_plants)~ca$map))
c_num_ft <- (lm(sqrt(ca$num_plants)~ca$canopy))
c_num_add <- (lm(sqrt(ca$num_plants)~ca$map+ca$canopy))
c_num_int <- (lm(sqrt(ca$num_plants)~ca$map*ca$canopy))


# Plot height as function of map
png(paste0("C:/Users/rache/Dropbox/Doctoral Projects/Data/Figures/",pftname,"_scatterplots.png"), 
    width = 7, height = 9, units = "in", res = 300)
par(mfrow = c(3,2))
par(mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(sqrt(avg_height)~map, data = ms, pch = 16,
     col = c(rgb(0,0.5,0,0.5),rgb(0.6,0.3,0.2,0.5))[as.factor(microsite)],
     main = paste0(pftname,"~Microsite: height"), xlab = "Mean annual precipitation (mm)",
     ylab = "sqrt(Average height)")
pickmodel(list(m_ht_map, m_ht_ft,  m_ht_add, m_ht_int), x)
plot(sqrt(avg_height)~map, data = ca, pch = 16,
     col = c(rgb(0,0.5,0,0.5),rgb(0.6,0.3,0.2,0.5))[as.factor(canopy)],
     main = paste0(pftname,"~Canopy cover: height"), xlab = "Mean annual precipitation (mm)",
     ylab = "sqrt(Average height)")
pickmodel(list(c_ht_map, c_ht_ft,  c_ht_add, c_ht_int), x)

# Plot number of plants as function of map
par(mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(sqrt(num_plants)~map, data = ms, pch = 16,
     col = c(rgb(0,0.5,0,0.5),rgb(0.6,0.3,0.2,0.5))[as.factor(microsite)],
     main = paste0(pftname,"~Microsite: number"), xlab = "Mean annual precipitation (mm)",
     ylab = "sqrt(Number of plants)")
pickmodel(list(m_num_map, m_num_ft,  m_num_add, m_num_int), x)
plot(sqrt(num_plants)~map, data = ca, pch = 16,
     col = c(rgb(0,0.5,0,0.5),rgb(0.6,0.3,0.2,0.5))[as.factor(canopy)],
     main = paste0(pftname,"~Canopy cover: number"), xlab = "Mean annual precipitation (mm)",
     ylab = "sqrt(Number of plants)")
pickmodel(list(c_num_map, c_num_ft,  c_num_add, c_num_int), x)

# Plot cover as function of map
par(mar = c(2,2,2,1), mgp = c(1,0.1,0), tcl = 0.1)
plot(sqrt(cover)~map, data = ms, pch = 16,
     col = c(rgb(0,0.5,0,0.5),rgb(0.6,0.3,0.2,0.5))[as.factor(microsite)],
     main = paste0(pftname,"~Microsite: cover"), xlab = "Mean annual precipitation (mm)",
     ylab = "sqrt(cover)")
pickmodel(list(m_cov_map, m_cov_ft,  m_cov_add, m_cov_int), x)
plot(sqrt(cover)~map, data = ca, pch = 16,
     col = c(rgb(0,0.5,0,0.5),rgb(0.6,0.3,0.2,0.5))[as.factor(canopy)],
     main = paste0(pftname,"~Canopy cover: cover"), xlab = "Mean annual precipitation (mm)",
     ylab = "sqrt(cover)")
pickmodel(list(c_cov_map, c_cov_ft,  c_cov_add, c_cov_int), x)
dev.off()

}

# Run plot_results function

x <- seq(167, 424, by = 1)

# All plants:
plot_results(pfts = c("CBG","CRG","CAG","PF","AF","SSH","WRG","WBG","SE","SUC"),x, pftname = "allplants")

# Pherbs
plot_results(pfts = c("CBG","CRG","WRG","WBG","PF","SE"), x, pftname = "Pherbs")

# Perennial grasses
plot_results(pfts=c("CBG","CRG","WRG","WBG"), x, pftname = "PG")

# Perennial bunch grasses
plot_results(pfts=c("CBG","WBG"), x, pftname = "PBG")

# Perennial rhizomatous grasses
plot_results(pfts=c("CRG","WRG"), x, pftname = "PRG")

# C4 Pgrasses
plot_results(pfts=c("WRG","WBG"), x, pftname = "C4Pgrass")

# C3 Pgrasses
plot_results(pfts=c("CRG","CBG"), x, pftname = "C3Pgrass")

# Perennial forbs
plot_results(pfts=c("PF"), x, pftname = "PF")

# Annual forbs
plot_results(pfts=c("AF"), x, pftname = "AF")

# CAG
plot_results(pfts=c("CAG"), x, pftname = "CAG")

# Bring in target shrub heights
target <- read.csv("C:/Users/rache/Dropbox/Doctoral Projects/Data/Final Data 2022/targetshrubs.csv")

# Restrict to target shrubs with actual neighbors
for (i in 1:nrow(nogaps)){
  target <- target[-which(target$plot == nogaps$plot[i] & 
                        target$transect == nogaps$transect[i] &
                        target$point == nogaps$point[i]),]
}

# Merge with plots
target1 <- merge(plots, target)
target1$canopy_size <- target1$length/2*target1$width/2*pi

plot(target1$height~target1$map, main = "Target shrub height by precipitation",
     xlab = "Mean annual precipitation (mm)", ylab = 'Height (cm)')
plot(target1$canopy_size~target1$map, main = "Target shrub canopy by precipitation",
     xlab = "Mean annual precipitation (mm)", ylab = 'Canopy area (cm2)')


##################################
# Now create relative-interaction index for each condition, summarized by plot for analysis of precipitation
# 
# # Convert NA's to 0 for rii
# rii <- community
# rii[is.na(rii$height),]$height <- 0
# rii[is.na(rii$number),]$number <- 0
# 
# # First summarize by microsite
# rii_microsite <- rii %>%
#   group_by(plot, transect, point, microsite) %>%
#   summarize(num_plants = sum(number, na.rm = T),
#             avg_height = mean(height, na.rm = T)) %>%
#   pivot_wider(names_from = microsite, values_from = c(num_plants, avg_height)) %>%
#   mutate(rii_num = (num_plants_C - num_plants_I)/(num_plants_C + num_plants_I), 
#          rii_height = (avg_height_C - avg_height_I)/(avg_height_C + avg_height_I))
# # Convert NA's to 0 in rii indexes
# rii_microsite[is.na(rii_microsite$rii_num),]$rii_num <- 0
# rii_microsite[is.na(rii_microsite$rii_height),]$rii_height <- 0
# 
# # Then summarize by canopy condition
# rii_canopy <- rii %>%
#   group_by(plot, transect, point, canopy) %>%
#   summarize(num_plants = sum(number, na.rm = T),
#             avg_height = mean(height, na.rm = T)) %>%
#   pivot_wider(names_from = canopy, values_from = c(num_plants, avg_height)) %>%
#   mutate(rii_num = (num_plants_Covered - num_plants_Open)/
#            (num_plants_Covered + num_plants_Open), 
#          rii_height = (avg_height_Covered - avg_height_Open)/
#            (avg_height_Covered + avg_height_Open))
# # Remove shrubs without both Open and Covered canopy conditions
# rii_canopy <- rii_canopy[complete.cases(rii_canopy[,c(1:7)]),]
# # Convert NA's to 0
# rii_canopy[is.na(rii_canopy$rii_num),]$rii_num <- 0
# rii_canopy[is.na(rii_canopy$rii_height),]$rii_height <- 0
# ```
# 
# 
# ```{r Test for normality of rii,message=FALSE, warning=FALSE, tidy=TRUE, include = FALSE, echo = FALSE}
# # Look at microsite
# par(mfrow = c(1,2), mar = c(2,2,2,1), tcl = 0.3, mgp = c(1,0,0))
# # Number of plants
# qqnorm(rii_microsite$rii_num)
# qqline(rii_microsite$rii_num)
# # Height
# qqnorm(rii_microsite$rii_height)
# qqline(rii_microsite$rii_height)
# 
# # Look at canopy condition
# par(mfrow = c(1,2), mar = c(2,2,2,1), tcl = 0.3, mgp = c(1,0,0))
# # Number of plants
# qqnorm(rii_canopy$rii_num)
# qqline(rii_canopy$rii_num)
# # Height
# qqnorm(rii_canopy$rii_height)
# qqline(rii_canopy$rii_height)
# ```
# 
# 
# ```{r Add precipitation data, message=FALSE, warning=FALSE, tidy=TRUE, include = FALSE, echo = FALSE}
# # Add precipitation data to each abundance data frame
# rii_microsite1 <- merge(rii_microsite, plots[,c(1,7)], 
#                         by.x = "plot",by.y ="ident")
# rii_canopy1 <- merge(rii_canopy, plots[,c(1,7)], 
#                      by.x = "plot",by.y ="ident")
# ```
# 
# 
