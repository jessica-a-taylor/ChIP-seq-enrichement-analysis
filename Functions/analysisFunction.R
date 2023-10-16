analysisFunction <- function(sampleGenes, nextflowOutput) {
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
    jobRunScript("Script for analysis.R", name= set, importEnv = TRUE)
  }
}