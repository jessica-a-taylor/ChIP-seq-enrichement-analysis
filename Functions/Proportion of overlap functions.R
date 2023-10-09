# Determine the proportion of each R-gene promotor overlapping with a significant peak.
promoterEnrichmentFunction <- function (allOverlaps, nextflowOutput, genomicData, promoterEnrichment, set) {
  
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
  promoterEnrichment[[set]] <- proportionDF
  return(promoterEnrichment)
}

# Determine the proportion of each R-gene body overlapping with a significant peak.
genebodyEnrichmentFunction <- function (allOverlaps, nextflowOutput, genomicData, genebodyEnrichment, set) {
  
  proportionDF <- data.frame(Gene = character(),
                             `Mod.TF` = character(),
                             Proportion = numeric())
  
  for (n in names(allOverlaps)) {
    for (mod in unique(nextflowOutput[, "Mod.TF"])) {
      peakWidths <- c()
      
      if (nrow(allOverlaps[[n]][[mod]]) >= 1) {
        for (row in 1:nrow(allOverlaps[[n]][[mod]])) {
          if (overlapsFunction(as.numeric(allOverlaps[[n]][[mod]][row,"start"]), as.numeric(allOverlaps[[n]][[mod]][row,"end"]),
                               as.numeric(geneRegions[["Gene20"]][geneRegions[["Gene20"]]$Gene==n,"start"]),
                               as.numeric(geneRegions[["Gene100"]][geneRegions[["Gene100"]]$Gene==n,"end"]))==TRUE) {
            
            peakWidths <- append(peakWidths, allOverlaps[[n]][[mod]][row,"width"])
          } else peakWidths <- peakWidths
        }
        
        if ((as.numeric(geneRegions[["Gene100"]][geneRegions[["Gene100"]]$Gene==n,"end"])-as.numeric(geneRegions[["Gene20"]][geneRegions[["Gene20"]]$Gene==n,"start"])) >= sum(peakWidths)) {
          proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                         `Mod.TF` = mod,
                                                         Proportion = sum(peakWidths)/(as.numeric(geneRegions[["Gene100"]][geneRegions[["Gene100"]]$Gene==n,"end"])-as.numeric(geneRegions[["Gene20"]][geneRegions[["Gene20"]]$Gene==n,"start"])),
                                                         ExpressionLevel = str_match(set, "^R-gene (.*)$")[,2])) 
          
        } else if ((as.numeric(geneRegions[["Gene100"]][geneRegions[["Gene100"]]$Gene==n,"end"])-as.numeric(geneRegions[["Gene20"]][geneRegions[["Gene20"]]$Gene==n,"start"]))  < sum(peakWidths)) {
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
  genebodyEnrichment[[set]] <- proportionDF
  return(genebodyEnrichment)
}


# Determine the proportion of each gene region overlapping with a significant peak.
proportionFunction <- function (geneRegions, allPeaks, nextflowOutput) {
  
  proportion <- hash()
  
  allModifications <- unique(nextflowOutput$Mod.TF)
  allRegions <- names(geneRegions)

  for (gene in names(allPeaks)) {
    for (region in allRegions) {
      for (mod in allModifications) {

        proportion[[gene]][[region]][[mod]] <- data.frame()
        
        if (nrow(allPeaks[[gene]][[region]][[mod]]) >= 1) {
          peakOverlaps <- c()
          
          for (row in 1:nrow(allPeaks[[gene]][[region]][[mod]])) {
            peakOverlaps <- append(peakOverlaps, newOverlapsFunction(as.numeric(allPeaks[[gene]][[region]][[mod]][row, "start"]), 
                                                                     as.numeric(allPeaks[[gene]][[region]][[mod]][row, "end"]),
                                                                     as.numeric(geneRegions[[region]][[gene]]$start),
                                                                     as.numeric(geneRegions[[region]][[gene]]$end)))
          }
          proportion[[gene]][[region]][[mod]] <- rbind(proportion[[gene]][[region]][[mod]], 
                                                       data.frame(Gene = gene,
                                                                  Region = region,
                                                                  `Mod.TF` = mod,
                                                                  Proportion = sum(peakOverlaps)/(as.numeric(geneRegions[[region]][[gene]]$width))))
        }
        else proportion[[gene]][[region]][[mod]] <- proportion[[gene]][[region]][[mod]]
      }
    }
  }
  return(proportion)
}