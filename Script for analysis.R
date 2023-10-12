library(readxl)
library(karyoploteR)
library(rtracklayer)
library(dplyr)
library(stringr)
library(hash)
library(sets)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(data.table)
library(grid)
library(readr)
library(rstudioapi)
library(ggplot2)
library(ggpubr)

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
  proportion <- geneRegionAxisLocations(proportion, geneRegions)
  
  # Store final results in 'proportionOfOverlap'.
  proportionOfOverlap[[set]] <- proportion
  
  print(set)
}
  
# Merge all data on the proportion of overlap per gene region into a single dataframe.
allProportions <- data.frame()

for (set in names(proportionOfOverlap)) {
  for (gene in names(proportionOfOverlap[[set]])) {
    for (region in names(proportionOfOverlap[[set]][[gene]])) {
      for (mod in names(proportionOfOverlap[[set]][[gene]][[region]])) {
        allProportions <- rbind(allProportions, proportionOfOverlap[[set]][[gene]][[region]][[mod]])
      }
    }
  }
}

#write.csv(geneFrequency, paste("PlantExp data\\allFrequencies.csv", sep = "")) 
write.csv(allProportions, paste("PlantExp data\\allProportions.csv", sep = "")) 