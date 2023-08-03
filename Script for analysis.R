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
source("Functions\\Proportion of overlap functions.R")
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")
source("Functions\\AxisGroup column.R")

# Create a hash for storing the proportion of each gene region overlapping with a significant peak.
sampleGenesProportionsPerRegion <- hash()

# Generate a list of the number of genes in each set.
geneCount <- data.frame()

for (set in names(sampleGenes)) {
    
  geneSet <- sampleGenes[[set]]

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
  
  # Determine the proportion of each gene overlapping with a significant peak.
  proportionPerGene <- hash()
  
  if (set %in% c("R-gene Low Expression","R-gene No Expression"  )) {
    proportionPerGene <- proportionPerGeneFunction(allOverlaps, nextflowOutput, genomicData, proportionPerGene, set)
  }
  
  # Determine the proportion of each gene region overlapping with a significant peak.
  geneRegions <- getGeneCoordinates(geneSet)
  
  proportionPerRegion <- proportionPerRegionFunction(geneRegions, allOverlaps, nextflowOutput)
  
  # Add a column to 'proportionPerRegion' with the numbers for 
  # each gene region that will correspond with their position on the x axis.
  proportionPerRegion <- geneRegionAxisLocations(proportionPerRegion, geneRegions)
  
  # Add a column to 'proportionPerRegion' identifying the gene set.
  proportionPerRegion <- expressionColumn(proportionPerRegion, set)
  
  # Store final results in 'sampleGenesProportionsPerRegion'.
  sampleGenesProportionsPerRegion[[set]] <- proportionPerRegion
  
  print(set)
}

# Merge all data on the proportion of overlap per gene into a single dataframe.
allProportionsPerGene <- data.frame()

for (set in names(proportionPerGene)) {
  allProportionsPerGene <- rbind(allProportionsPerGene, proportionPerGene[[set]])
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
  write.csv(allProportionsPerGene, paste("PlantExp data\\Normalised\\allProportionsPerGene.csv", sep = "")) 
  write.csv(allProportionsPerRegion, paste("PlantExp data\\Normalised\\allProportionsPerRegion.csv", sep = "")) 
} else if (normalised == FALSE) {
  write.csv(geneFrequency, paste("PlantExp data\\Non-normalised\\allFrequencies.csv", sep = "")) 
  write.csv(allProportionsPerGene, paste("PlantExp data\\Non-normalised\\allProportionsPerGene.csv", sep = "")) 
  write.csv(allProportionsPerRegion, paste("PlantExp data\\Non-normalised\\allProportionsPerRegion.csv", sep = "")) 
}

