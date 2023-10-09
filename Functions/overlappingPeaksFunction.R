overlappingPeaksFunction <- function(geneRegions, allModifications, allRegions, peaksPerGene, allPeaks, gene) {

  # Initialise an empty hash of dataframes.
  overlappingPeaks <- hash()

  for (mod in allModifications) {

    # Filter for each modification/TF.
    peaksPerModification <- peaksPerGene[peaksPerGene$Mod.TF == mod,]

    for (region in names(geneRegions)) {
      overlappingPeaks[[gene]][[region]][[mod]] <- data.frame()
      allPeaks[[gene]][[region]][[mod]] <- data.frame()
    }
    
    # Populate the hash
    # For each row in peaksPerModification,
    for (row in 1:nrow(peaksPerModification)) {
      
      # Get the start and end coordinates of the gene region,
      for (region in allRegions) {
        regionStart <- geneRegions[[region]][[gene]][, "start"]
        regionEnd <- geneRegions[[region]][[gene]][, "end"]
        
        # Get the peaks that overlap with each gene region.
        if (gene %in% names(geneRegions[[region]]) & nrow(peaksPerModification) > 0) {
          if (overlapsFunction(peaksPerModification[row, "start"], peaksPerModification[row, "end"],
                               regionStart, regionEnd) == TRUE) {
            
            overlappingPeaks[[gene]][[region]][[mod]] <- rbind(overlappingPeaks[[gene]][[region]][[mod]], peaksPerModification[row,])
            
          }
        }
      }
    }
    # Merge overlapping peaks for the same modification/TF.
    for (region in allRegions) {
      # Generate overlapSets as a list of single-item sets
      # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
      
      overlapSets <- list()
      
      numberOfOverlappingPeaks <- nrow(overlappingPeaks[[gene]][[region]][[mod]])
      
      if (numberOfOverlappingPeaks > 0) {
        for (row in 1:numberOfOverlappingPeaks) {
          
          overlapSets[[row]] <- list()
          
          overlapSets[[row]][["set"]] <- sets::set(as.numeric(row))
          overlapSets[[row]][["start"]] <- overlappingPeaks[[gene]][[region]][[mod]][row,"start"]
          overlapSets[[row]][["end"]] <- overlappingPeaks[[gene]][[region]][[mod]][row,"end"]
        }
        
        k <- 1
        
        # For each gene co-ordinate comparison [k, l]
        while (k < length(overlapSets)) {

          foundAnOverlap <- TRUE
          
          while (foundAnOverlap == TRUE) {
            
            foundAnOverlap <- FALSE

            l <- k + 1
            
            while (l <= length(overlapSets)) {

              # If the co-ordinate ranges overlap
              if (overlapsFunction(overlapSets[[k]][["start"]],overlapSets[[k]][["end"]],
                                   overlapSets[[l]][["start"]],overlapSets[[l]][["end"]]) == TRUE) {
                
                # Merge the two sets, replacing the old ones
                overlapSets[[k]][["start"]] <- min(overlapSets[[k]][["start"]], overlapSets[[l]][["start"]])
                overlapSets[[k]][["end"]] <- max(overlapSets[[k]][["end"]], overlapSets[[l]][["end"]])
                overlapSets[[k]][["set"]] <- set_union(overlapSets[[k]][["set"]], overlapSets[[l]][["set"]])
                
                overlapSets <- overlapSets[-c(l)]
              }
              else {
                l <- l + 1
              }
            }
          }
          k <- k + 1
        }
        # Find the range of the merged peaks.
        for (l in length(overlapSets)) {
          
          start <- overlapSets[[l]][["start"]]
          end <- overlapSets[[l]][["end"]]
          
          allPeaks[[gene]][[region]][[mod]] <- rbind(allPeaks[[gene]][[region]][[mod]], 
                                                     data.frame(gene = gene,
                                                                seqnames = overlappingPeaks[[gene]][[region]][[mod]]$seqnames[1],
                                                                start = start,
                                                                end = end,
                                                                width = end - start,
                                                                ranges = paste(start, "-", end, sep = ""),
                                                                region = region,
                                                                `Mod.TF` = mod))
        }
      }
      else next
    }
  }
  return(allPeaks)
}