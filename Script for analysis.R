source("Functions\\Coordinates per gene region.R")
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Overlaps functions.R")
source("Functions\\peakOverlaps.R")
source("Functions\\proportionFunctions.R")
source("Functions\\AxisGroup column.R")

genomicData <- as.data.frame(read_csv("Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Create a hash for storing the proportion of each gene region overlapping with a significant peak.
proportionOfOverlap <- hash()

# Generate a list of the number of genes in each set.
geneCount <- data.frame()

for (set in names(sampleGenes)) {
  
  wb <- createWorkbook()
    
  geneSet <- sampleGenes[[set]]
  
  # Get the coordinates for each genomic region.
  geneRegions <- getGeneCoordinates(geneSet, genomicData)

  # For each modification/TF, get the peaks overlapping each gene region.
  allPeaks <- peakOverlaps(geneSet, geneRegions, nextflowOutput, set)
  
  geneCount <- rbind(geneCount, data.frame(GeneSet = set,
                                           GeneCount = length(names(allPeaks))))

  
  # Determine the proportion of each gene region overlapping with a significant peak.
  proportion <- proportionFunction(geneRegions, allPeaks, nextflowOutput)
  
  # Add a column to 'proportion' with the numbers for 
  # each gene region that will correspond with their position on the x axis.
  wb <- geneRegionAxisLocations(proportion, geneRegions, wb)
  
  # Store final results in 'proportionOfOverlap'.
  #proportionOfOverlap[[set]] <- proportion
  
  print(set)
  
  saveWorkbook(wb, paste(set, ".xlsx", sep = ""), overwrite = TRUE)
}