library(readxl)
library(karyoploteR)
library(rtracklayer)
library(dplyr)
library(stringr)
library(hash)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(data.table)
library(grid)
library(readr)
library(rstudioapi)
library(ggplot2)
library(ggpubr)

source("Functions\\Overlaps functions.R")
source("Functions\\Peaks per gene.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Proportion of gene region.R")
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")
source("Functions\\AxisGroup column.R")

# Create a hash for storing the proportion of each gene region overlapping with a significant peak.
sampleGenesPlantExpProportions <- hash()

# Generate a list of the number of genes in each set.
geneCount <- data.frame()

for (set in names(sampleGenesPlantExp)) {
    
  geneSet <- sampleGenesPlantExp[[set]]

  # Create a hash containing the significant peaks in each gene. 
  allPeaks <- PeaksPerGene(geneSet, nextflowOutput)
  
  # Create a hash containing hashes with the coordinates of significant peaks from each experiment for each gene.
  genePeaks <- peakOccurrences(allPeaks, nextflowOutput)
  
  rm(allPeaks)
  
  # For each gene in the current set of genes, merge overlapping peaks.
  allOverlaps <- mergeOverlappingPeaks(genePeaks, nextflowOutput)
  
  geneCount <- rbind(geneCount, data.frame(GeneSet = set,
                                           GeneCount = length(names(allOverlaps))))
  rm(genePeaks)
  
  # Determine the proportion of each gene region overlapping with a significant peak.
  geneRegions <- getGeneCoordinates(geneSet)
  
  proportionPerRegion <- proportionsFunction(geneRegions, allOverlaps, nextflowOutput)
  
  # Add a column to 'proportionPerRegion' with the numbers for 
  # each gene region that will correspond with their position on the x axis.
  proportionPerRegion <- geneRegionAxisLocations(proportionPerRegion, geneRegions)
  
  # Add a column to 'proportionPerRegion' identifying the gene set.
  proportionPerRegion <- expressionColumn(proportionPerRegion, set)
  
  # Store final results in 'sampleGenesPlantExpProportions'.
  sampleGenesPlantExpProportions[[set]] <- proportionPerRegion
  
  print(set)
}
  
# Merge all data from all sample gene sets into one big dataframe.
allResultsProportions <- data.frame()

for (set in names(sampleGenesPlantExpProportions)) {
  df <- sampleGenesPlantExpProportions[[set]]
  allResultsProportions <- rbind(allResultsProportions, df)
}

# Generate a list of the number of genes modified in the control- and R-gene set &
# determine the variance in the enrichment for each modification/TF in the control- and R-gene set.

source("Functions\\% enriched genes.R")

geneFrequency <- data.frame(Region = rep(unique(allResultsProportions$Region), times = 4*length(unique(allResultsProportions$Mod.TF))),
                            Mod.TF = rep(unique(allResultsProportions$Mod.TF), each = 4*length(unique(allResultsProportions$Region))),
                            GeneSet = rep(unique(allResultsProportions$GeneSet), each = length(unique(allResultsProportions$Region))),
                            Count = rep(c(0,0), times = 2*length(unique(allResultsProportions$Region))),
                            Enrichment.mean = rep(0, times = 4*length(unique(allResultsProportions$Region))),
                            Enrichment.variance = rep(0, times = 4*length(unique(allResultsProportions$Region))))

geneFrequency <- frequenciesFunction(allResultsProportions, geneFrequency, geneCount)

# Add a column to 'proportionPerRegion' with the numbers for each gene region that will correspond with their position on the x axis.
geneFrequency <- geneRegionAxisLocations(geneFrequency, geneRegions)

if (normalised == TRUE) {
  write.csv(geneFrequency, paste("PlantExp data\\Normalised\\allResultsFrequencies.csv", sep = "")) 
  write.csv(allResultsProportions, paste("PlantExp data\\Normalised\\allResultsProportions.csv", sep = "")) 
} else if (normalised == FALSE) {
  write.csv(geneFrequency, paste("PlantExp data\\Non-normalised\\allResultsFrequencies.csv", sep = "")) 
  write.csv(allResultsProportions, paste("PlantExp data\\Non-normalised\\allResultsProportions.csv", sep = "")) 
}

