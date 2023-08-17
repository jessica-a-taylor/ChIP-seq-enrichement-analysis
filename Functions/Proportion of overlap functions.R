source("Functions\\Overlaps functions.R")

# Determine the proportion of each gene overlapping with a significant peak.
proportionPerGeneFunction <- function (allOverlaps, nextflowOutput, genomicData, proportionPerGene, set) {
  
  proportionDF <- data.frame(Gene = character(),
                             `Mod.TF` = character(),
                             Proportion = numeric())
  
  for (n in names(allOverlaps)) {
    for (mod in unique(nextflowOutput[, "Mod.TF"])) {
      peakWidths <- c()
      
      if (nrow(allOverlaps[[n]][[mod]]) >= 1) {
        for (row in 1:nrow(allOverlaps[[n]][[mod]])) {
          if (overlapsFunction(as.numeric(allOverlaps[[n]][[mod]][row,"start"]), as.numeric(allOverlaps[[n]][[mod]][row,"end"]),
                               geneRegions[["Promotor1000"]][geneRegions[["Promotor1000"]]$Gene==n,"start"],
                               geneRegions[["Promotor1000"]][geneRegions[["Promotor1000"]]$Gene==n,"end"])==TRUE) {
            
            peakWidths <- append(peakWidths, allOverlaps[[n]][[mod]][row,"width"])
          } else peakWidths <- peakWidths
        }
        
        if (1000 >= sum(peakWidths)) {
          proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                         `Mod.TF` = mod,
                                                         Proportion = sum(peakWidths)/1000,
                                                         ExpressionLevel = str_match(set, "^R-gene (.*)$")[,2])) 
          
        } else if (1000 < sum(sum(peakWidths))) {
          proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                         `Mod.TF` = mod,
                                                         Proportion = 1,
                                                         ExpressionLevel = str_match(set, "^R-gene (.*)$")[,2]))
        }
      }
      else if (nrow(allOverlaps[[n]][[mod]]) < 1) {
        proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                       `Mod.TF` = mod,
                                                       Proportion = 0,
                                                       ExpressionLevel = str_match(set, "^R-gene (.*)$")[,2]))
      }
    }
  } 
  proportionPerGene[[set]] <- proportionDF
  return(proportionPerGene)
}


# Determine the proportion of each gene region overlapping with a significant peak.
proportionPerRegionFunction <- function (geneRegions, allOverlaps, nextflowOutput) {

  proportionPerRegion <- hash()
  
  for (r in names(geneRegions)) {
    proportionDF <- data.frame(Gene = character(),
                               Region= character(), 
                               `Mod.TF` = character(),
                               Proportion = numeric())
    
    if (length(names(allOverlaps)) >= 1) {
      for (mod in unique(nextflowOutput[, "Mod.TF"])) {
        
        for (n in names(allOverlaps)) {
          
          peakOverlaps <- c()
          
          if (nrow(allOverlaps[[n]][[mod]]) >= 1 & n %in% geneRegions[[r]]$Gene == TRUE) {
            
            for (row in 1:nrow(allOverlaps[[n]][[mod]])) {
              peakOverlaps <- append(peakOverlaps, newOverlapsFunction(as.numeric(allOverlaps[[n]][[mod]][row, "start"]), as.numeric(allOverlaps[[n]][[mod]][row, "end"]),
                                                                     as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$start), as.numeric(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$end)))
            }
            proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                             Region = r,
                                                             `Mod.TF` = mod,
                                                             Proportion = sum(peakOverlaps)/(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width)))
          }
          else proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                              Region = r,
                                                              `Mod.TF` = mod,
                                                              Proportion = 0))
        }
      }
    } else proportionDF <- proportionDF
    
    proportionPerRegion[[r]] <- proportionDF
  }
  
  # Collect all hashes into a single dataframe.
  DF <- data.frame()
  
  for (r in names(geneRegions)) {
    DF <- rbind(DF, proportionPerRegion[[r]])
  }
  return(DF)
}