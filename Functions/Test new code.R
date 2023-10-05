library(hash)
library(sets)

Rprof(line.profiling = TRUE)

overlappingPeaks <- hash()
mergedPeaks <- data.frame()

filteringTime <- 0
populatingTime <- 0
mergeTime <- 0

# Filter for each modification.
# nextflowOutput <- nextflowOutput[nextflowOutput$Mod.TF == mod,]

# Running unique() on a dataframe will require iterating through the whole thing,
# so do it once instead
allModifications <- unique(nextflowOutput$Mod.TF)
allRegions <- names(geneRegions)

executionStart <- Sys.time()

i <- 0

totalIterations <- length(geneSet$Gene)

for (gene in geneSet$Gene) {
  i = i + 1
  print(paste(gene, " (", i, "/", totalIterations, ")"))
  
  # Filter for each gene.
  df <- geneSet[geneSet$Gene == gene,]
  
  filterStart <- Sys.time()
  
  # Filter for peaks overlapping each gene, promotor and surrounding intergenic regions.
  peaks_df <- nextflowOutput[which(nextflowOutput[, "start"] > df[1, "start"] - 5000 &
                                   nextflowOutput[, "end"] < df[1, "end"] + 5000 &
                                   as.numeric(nextflowOutput[, "seqnames"]) == as.numeric(df[1, "seqnames"])),]
  
  filterEnd <- Sys.time()
  
  filteringTime <- filteringTime + filterEnd - filterStart
  
  dataLength <- nrow(peaks_df)
  
  allRows <- 1:dataLength
  
  # Initialise an empty hash of dataframes
  for (mod in allModifications) {
    for (region in names(geneRegions)) {
      overlappingPeaks[[gene]][[region]][[mod]] <- data.frame()
    }
  }
  
  populateStart <- Sys.time()
  
  # Populate the hash
  for (row in allRows) {
    for (region in allRegions) {
      regionStart <- geneRegions[[region]][[gene]][, "start"]
      regionEnd <- geneRegions[[region]][[gene]][, "end"]
      
      if (gene %in% names(geneRegions[[region]]) & dataLength > 0) {
        if (overlapsFunction(peaks_df[row, "start"],
                             peaks_df[row, "end"],
                             regionStart,
                             regionEnd) == TRUE) {
          
          overlappingPeaks[[gene]][[region]][[mod]] <- rbind(overlappingPeaks[[gene]][[region]][[mod]], peaks_df[row,])
          
        }
      }
    }
  }
  
  populateEnd <- Sys.time()
  populatingTime <- populatingTime + populateEnd - populateStart
  
  mergeStart <- Sys.time()
  for (mod in allModifications) {
    # Merge overlapping peaks for the same modification/TF.
    
    for (region in names(overlappingPeaks[[gene]])) {
      # Generate overlapSets as a list of single-item sets
      # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
      
      overlapSets <- list()
      
      overlappingPeaksLength <- nrow(overlappingPeaks[[gene]][[region]][[mod]])
      
      if (overlappingPeaksLength > 0) {
        for (row in 1:overlappingPeaksLength) {
          
          overlapSets[[row]] <- list()
          
          overlapSets[[row]][["set"]] <- sets::set(as.numeric(row))
          overlapSets[[row]][["start"]] <- overlappingPeaks[[gene]][[region]][[mod]][row,"start"]
          overlapSets[[row]][["end"]] <- overlappingPeaks[[gene]][[region]][[mod]][row,"end"]
        }
        
        k <- 1
        
        # For each gene co-ordinate comparison [k, l]
        while (k <= length(overlapSets)) {
          foundAnOverlap <- TRUE
          
          while (foundAnOverlap == TRUE) {
            
            foundAnOverlap <- FALSE
            
            l <- k + 1
            
            while (l <= length(overlapSets)) {
              
              # If the co-ordinate ranges overlap
              if (overlapsFunction(
                overlapSets[[k]][["start"]],
                overlapSets[[k]][["end"]],
                overlapSets[[l]][["start"]],
                overlapSets[[l]][["end"]]
              ) == TRUE) {
                
                foundAnOverlap <- TRUE
                  
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
      }
      else  next
      
      # Find the range of the merged peaks.
      for (l in length(overlapSets)) {
        print(paste(mod, region, l))
        
        start <- overlapSets[[l]][["start"]]
        end <- overlapSets[[l]][["end"]]
        
        mergedPeaks <- rbind(
          mergedPeaks,
          data.frame(
            gene = gene,
            seqnames = overlappingPeaks[[gene]][[region]][[mod]]$seqnames[1],
            start = start,
            end = end,
            width = end - start,
            ranges = paste(start, "-", end, sep = ""),
            region = region,
            `Mod.TF` = mod
          )
        )
      }
      
    }
  }
  mergeEnd <- Sys.time()
  
  mergeTime <- mergeTime + mergeEnd - mergeStart
}

executionEnd <- Sys.time()

print(paste(
  "Completed in ",
  as.numeric(executionEnd - executionStart,units="secs"),
  ". Time spent filtering: ",
  as.numeric(filteringTime,units="secs"),
  "; time spent populating: ",
  as.numeric(populatingTime,units="secs"),
  "; time spent merging: ",
  as.numeric(mergeTime,units="secs")
))

write.csv(mergedPeaks, file = "mergedPeaks.csv")