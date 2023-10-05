library(hash)
library(sets)
library(stringr)

allPeaks <- hash()

if (nrow(geneSet) >= 1) {
  for (row in 1:nrow(geneSet)) {
    
    # Select rows that are within the range of each gene and on the same chromosome.
    selectedRows <- c(which(nextflowOutput[,"start"] > geneSet[row, "start"]-5000 & nextflowOutput[,"end"] < geneSet[row, "end"]+5000 & nextflowOutput[,"seqnames"] == as.numeric(geneSet[row, "seqnames"])))
    allPeaks[[geneSet[row,"Gene"]]] <- nextflowOutput[selectedRows,]
  }
}
return(allPeaks)


# Function for creating a hash containing hashes with the coordinates of significant peaks from each Mod.TF for each gene.
genePeaks <- hash()

for (n in names(allPeaks)) {
  peakHash <- hash()
  
  for (mod in unique(nextflowOutput[, "Mod.TF"])) {
    peakHash[[mod]] <- allPeaks[[n]][which(allPeaks[[n]][, "Mod.TF"]==mod),]
  }
  genePeaks[[n]] <- peakHash
}
return(genePeaks)


# Function for creating a hash with the overlapping occurrences of each modification merged together.
allOverlaps <- hash()

# For each Mod.TF...
for (n in names(genePeaks)) {
  peakOverlaps <- hash()
  
  for (mod in unique(nextflowOutput[, "Mod.TF"])) {
    
    # Generate overlapSets as a list of single-item sets
    # eg, [ {1}, {2}, {3}, {4}, {5}, {6} ]
    overlapSets <- list()
    if (nrow(genePeaks[[n]][[mod]])>0) {
      
      for (r in 1:nrow(genePeaks[[n]][[mod]])) {
        overlapSets <- append(overlapSets, list(sets::set(as.numeric(r))))
      }
      # For each gene co-ordinate comparison [k, l]
      for (k in 1:nrow(genePeaks[[n]][[mod]])) {
        for (l in 1:k) {
          
          # If the co-ordinate ranges overlap
          if (overlapsFunction(genePeaks[[n]][[mod]][k, "start"], genePeaks[[n]][[mod]][k, "end"], 
                               genePeaks[[n]][[mod]][l, "start"], genePeaks[[n]][[mod]][l, "end"])==TRUE) {
            
            # Find the indexes of the sets containing each range
            kIndex <- findItem(k, overlapSets)
            lIndex <- findItem(l, overlapSets)
            
            # No need to merge if the co-ordinate ranges are already in the same sets
            if (kIndex!=lIndex) {
              
              # If they are in different sets, merge the two sets, replacing the old ones
              newSet <- set_union(overlapSets[[kIndex]], overlapSets[[lIndex]])
              overlapSets <- overlapSets[-c(kIndex, lIndex)]
              overlapSets <- append(overlapSets, list(newSet))
            }
          }
        }
      }
    } 
    else next
    peakOverlaps[[mod]] <- overlapSets
  }
  allOverlaps[[n]] <- peakOverlaps
}


# Find the maximum range for the overlapping peaks.
for (n in names(allOverlaps)) {
  for (mod in unique(nextflowOutput[, "Mod.TF"])) {
    if (length(allOverlaps[[n]][[mod]])>0) {
      
      for (l in 1:length(allOverlaps[[n]][[mod]])) {
        peakStart <- c()
        peakEnd <- c() 
        
        for (o in allOverlaps[[n]][[mod]][l]) {
          peakStart <- append(peakStart, genePeaks[[n]][[mod]][as.numeric(o), "start"])
          peakEnd <- append(peakEnd, genePeaks[[n]][[mod]][as.numeric(o), "end"])
          
          allOverlaps[[n]][[mod]][l] <- paste(min(peakStart), max(peakEnd), sep = "-")
        }
      }
    }
  }
}
  
  
# Create dataframes with the information needed in the bed file.
for (n in names(genePeaks)) {
  for (mod in unique(nextflowOutput[, "Mod.TF"])) {
    df <- data.frame(seqnames = numeric(),
                     start = numeric(),
                     end = numeric(),
                     width = numeric(),
                     ranges = character(),
                     `Mod.TF` = character())
    
    if (length(allOverlaps[[n]][[mod]])>0) {
      for (l in 1:length(allOverlaps[[n]][[mod]])) {
        df <- rbind(df, data.frame(seqnames = genePeaks[[n]][[mod]][1,"seqnames"],
                                   start = str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,2],
                                   end = str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,3],
                                   width = as.numeric(str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,3]) - as.numeric(str_match(allOverlaps[[n]][[mod]][[l]], "^([0-9]+)-([0-9]+)$")[,2]),
                                   ranges = allOverlaps[[n]][[mod]][[l]],
                                   `Mod.TF` = mod))
      }
    }
    allOverlaps[[n]][[mod]] <- df
  }
}

print(allOverlaps[[n]][["H3K4me3"]])