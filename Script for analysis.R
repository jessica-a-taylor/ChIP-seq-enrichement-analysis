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

source("Functions\\Overlaps functions.R")
source("Functions\\peakOverlaps.R")
source("Functions\\Coordinates per gene region.R")
source("Functions\\Proportion of overlap functions.R")
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")
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
  allPeaks <- peakOverlaps(geneSet, geneRegions, nextflowOutput)
  
  geneCount <- rbind(geneCount, data.frame(GeneSet = set,
                                           GeneCount = length(names(allPeaks))))

  
  # Determine the proportion of each gene region overlapping with a significant peak.
  proportion <- proportionFunction(geneRegions, allPeaks, nextflowOutput)
  
  # Add a column to 'proportion' with the numbers for 
  # each gene region that will correspond with their position on the x axis.
  proportion <- geneRegionAxisLocations(proportion, geneRegions)
  
  # Add a column to 'proportion' identifying the gene set.
  proportion <- expressionColumn(proportion, set)
  
  # Store final results in 'proportionOfOverlap'.
  proportionOfOverlap[[set]] <- proportion
  
  print(set)
}

# Merge all data on the proportion of overlap per gene into a single dataframe.
allPromoterEnrichment <- data.frame()
allGenebodyEnrichment <- data.frame()

for (set in names(proportionPerGene)) {
  allPromoterEnrichment <- rbind(allPromoterEnrichment, promoterEnrichment[[set]])
  allGenebodyEnrichment <- rbind(allGenebodyEnrichment, genebodyEnrichment[[set]])
}
  
# Merge all data on the proportion of overlap per gene region into a single dataframe.
allProportionsPerRegion <- data.frame()

for (set in names(sampleGenesProportionsPerRegion)) {
  allProportionsPerRegion <- rbind(allProportionsPerRegion, sampleGenesProportionsPerRegion[[set]])
}

# Generate a list of the number of genes modified in the control- and R-gene set &
# determine the variance in the enrichment for each modification/TF in the control- and R-gene set.

source("Functions\\% enriched genes.R")

geneFrequency <- data.frame(Region = rep(unique(allProportionsPerRegion$Region), times = 4*length(unique(allProportionsPerRegion$Mod.TF))),
                            Mod.TF = rep(unique(allProportionsPerRegion$Mod.TF), each = 4*length(unique(allProportionsPerRegion$Region))),
                            GeneSet = rep(unique(allProportionsPerRegion$GeneSet), each = length(unique(allProportionsPerRegion$Region))),
                            Count = rep(c(0,0), times = 2*length(unique(allProportionsPerRegion$Region))),
                            Enrichment.mean = rep(0, times = 4*length(unique(allProportionsPerRegion$Region))),
                            Enrichment.variance = rep(0, times = 4*length(unique(allProportionsPerRegion$Region))))

geneFrequency <- frequenciesFunction(allProportionsPerRegion, geneFrequency, geneCount)

# Add a column to 'proportionPerRegion' with the numbers for each gene region that will correspond with their position on the x axis.
geneFrequency <- geneRegionAxisLocations(geneFrequency, geneRegions)

if (normalised == TRUE) {
  write.csv(geneFrequency, paste("PlantExp data\\Normalised\\allFrequencies.csv", sep = "")) 
  write.csv(allPromoterEnrichment, paste("PlantExp data\\Normalised\\promoterEnrichment.csv", sep = ""))
  write.csv(allGenebodyEnrichment, paste("PlantExp data\\Normalised\\genebodyEnrichment.csv", sep = "")) 
  write.csv(allProportionsPerRegion, paste("PlantExp data\\Normalised\\allProportionsPerRegion.csv", sep = "")) 
} else if (normalised == FALSE) {
  write.csv(geneFrequency, paste("PlantExp data\\Non-normalised\\allFrequencies.csv", sep = "")) 
  write.csv(allPromoterEnrichment, paste("PlantExp data\\Non-normalised\\promoterEnrichment.csv", sep = "")) 
  write.csv(allGenebodyEnrichment, paste("PlantExp data\\Non-normalised\\genebodyEnrichment.csv", sep = "")) 
  write.csv(allProportionsPerRegion, paste("PlantExp data\\Non-normalised\\allProportionsPerRegion.csv", sep = "")) 
}