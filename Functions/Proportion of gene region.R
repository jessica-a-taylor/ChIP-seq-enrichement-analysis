source("Functions\\Overlaps functions.R")

region = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
           "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
           "Downstream", "DownstreamIntergenic")

# Determine the proportion of each gene region overlapping with a significant peak.
proportionsFunction <- function (geneRegions, allOverlaps, data) {

  proportionPerRegion <- hash()
  
  for (r in names(geneRegions)) {
    proportionDF <- data.frame(Gene = character(),
                               Region= character(), 
                               `Mod/TF` = character(),
                               Proportion = numeric())
    
    if (length(names(allOverlaps)) >= 1) {
      for (mod in unique(data[, "Mod/TF"])) {
        
        for (n in names(allOverlaps)) {
          
          peakOverlaps <- c()
          
          if (nrow(allOverlaps[[n]][[mod]]) >= 1 & n %in% geneRegions[[r]]$Gene == TRUE) {
            
            for (row in 1:nrow(allOverlaps[[n]][[mod]])) {
              peakOverlaps <- append(peakOverlaps, newOverlapsFunction(as.numeric(allOverlaps[[n]][[mod]][row, "start"]), as.numeric(allOverlaps[[n]][[mod]][row, "end"]),
                                                                     as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$start), as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$end)))
            }
            proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                           Region = r,
                                                           `Mod/TF` = mod,
                                                           Proportion = sum(peakOverlaps)/(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width)*((geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width)/mean(geneRegions[[r]]$width))))
          }
          else proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                              Region = r,
                                                              `Mod/TF` = mod,
                                                              Proportion = 0))
        }
      }
    } else proportionDF <- proportionDF
    
    proportionPerRegion[[r]] <- proportionDF
  }
  
  # Collect all hashes into a single dataframe.
  DF <- data.frame(Gene = character(),
                   Region = character(),
                   Feature = character(),
                   Measure = numeric())
  
  for (r in region) {
    DF <- rbind(DF, proportionPerRegion[[r]])
  }
  return(DF)
}


# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(dataToUse, geneRegions) {
  grouping <- c(seq(from = -60, to = -20, by = 20), seq(from = 20, to = 100, by = 20), seq(from = 140, to = 160, by = 20))
  axisGroup <- c()
  
  for (c in 1:length(names(geneRegions))) {
    axisGroup <- append(axisGroup, rep(grouping[c], times = nrow(dataToUse[dataToUse$Region == names(geneRegions)[c],])))
  }
  dataToUse <- cbind(dataToUse, axisGroup)
  
  return(dataToUse) 
}


# Function to add a column for modression level.
modressionColumn <- function(dataToUse, test) {
  if (nrow(dataToUse) >= 1) {
    dataToUse <- cbind(dataToUse, data.frame(modression = rep(unlist(str_split(test, "_"))[2], times = nrow(dataToUse))))
  }
  else dataToUse <- dataToUse
  return(dataToUse)
}