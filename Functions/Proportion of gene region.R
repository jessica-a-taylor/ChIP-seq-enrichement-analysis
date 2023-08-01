source("Functions\\Overlaps functions.R")

region = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
           "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
           "Downstream", "DownstreamIntergenic")

# Determine the proportion of each gene region overlapping with a significant peak.
proportionsFunction <- function (geneRegions, allOverlaps, nextflowOutput) {

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
            if (normalised == FALSE) {
              proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                             Region = r,
                                                             `Mod.TF` = mod,
                                                             Proportion = sum(peakOverlaps)/(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width)))
            }
            else if (normalised == TRUE) {
              proportionDF <- rbind(proportionDF, data.frame(Gene = n,
                                                             Region = r,
                                                             `Mod.TF` = mod,
                                                             Proportion = sum(peakOverlaps)/(geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width)*((geneRegions[[r]][geneRegions[[r]]$Gene==n,]$width)/mean(geneRegions[[r]]$width))))
            }
            
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