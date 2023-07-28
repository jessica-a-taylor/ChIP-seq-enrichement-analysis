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
source("Functions\\% enriched genes.R")
source("Functions\\Proportion of gene region.R")
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")


if (analysis == "PlantExp data") {
  dataToAnalyse <- sampleGenesPlantExp
} else if (analysis == "RNA-seq data") {
  dataToAnalyse <- sampleGenesRNAseq
}

# Create a hash for storing the percentage of genes that overlap with a significant peak.
dataToAnalyseFrequencies <- hash()

# Create a hash for storing the proportion of each gene region overlapping with a significant peak.
dataToAnalyseProportions <- hash()

# Generate a list of the number of genes in each set.
geneCount <- data.frame()

for (normalised in c(TRUE, FALSE)) {
  for (test in names(dataToAnalyse)[grepl("_", names(dataToAnalyse))]) {
    
    geneSet <- dataToAnalyse[[test]]
    
    # Create a hash containing the significant peaks in each gene. 
    allPeaks <- PeaksPerGene(geneSet, nextflowOutput)
    
    # Create a hash containing hashes with the coordinates of significant peaks from each experiment for each gene.
    genePeaks <- peakOccurrences(allPeaks, nextflowOutput)
    
    rm(allPeaks)
    
    # For each gene in the current set of genes, merge overlapping peaks.
    allOverlaps <- mergeOverlappingPeaks(genePeaks, nextflowOutput)
    
    geneCount <- rbind(geneCount, data.frame(GeneSet = test,
                                             GeneCount = length(names(allOverlaps))))
    rm(genePeaks)
    
    # Determine the percentage of genes that overlap with a significant peak.
    # Determine the proportion of each gene region overlapping with a significant peak.
    geneRegions <- getGeneCoordinates(geneSet)
    
    percentageEnrichedGenes <- frequenciesFunction(geneRegions, allOverlaps, nextflowOutput)
    proportionPerRegion <- proportionsFunction(geneRegions, allOverlaps, nextflowOutput)
    
    # Add a column to 'percentageEnrichedGenes' & 'proportionPerRegion' with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    percentageEnrichedGenes <- geneRegionAxisLocations(percentageEnrichedGenes, geneRegions)
    proportionPerRegion <- geneRegionAxisLocations(proportionPerRegion, geneRegions)
    
    rm(geneRegions)
    
    # Add a column to 'percentageEnrichedGenes' & 'proportionPerRegion' with the current expression level.
    percentageEnrichedGenes <- expressionColumn(percentageEnrichedGenes, test)
    proportionPerRegion <- expressionColumn(proportionPerRegion, test)
    
    # Store final results in 'dataToAnalyseFrequencies' & 'dataToAnalyseProportions'.
    dataToAnalyseFrequencies[[test]] <- percentageEnrichedGenes
    dataToAnalyseProportions[[test]] <- proportionPerRegion
    
    print(test)
  }
  
  # Merge all data from all sample gene sets into one big dataframe.
  allResultsFrequencies <- data.frame()
  allResultsProportions <- data.frame()
  
  for (test in names(dataToAnalyseProportions)) {
    
    df <- dataToAnalyseFrequencies[[test]]
    df <- cbind(df, data.frame(dataToAnalyse = rep(unlist(str_split(test, "_"))[1], times = nrow(df))))
    
    allResultsFrequencies <- rbind(allResultsFrequencies, df)
    
    df <- dataToAnalyseProportions[[test]]
    df <- cbind(df, data.frame(dataToAnalyse = rep(unlist(str_split(test, "_"))[1], times = nrow(df))))
    
    allResultsProportions <- rbind(allResultsProportions, df)
  }
  
  # Calculate the mean proportion of overlap and the variance and add to 'allResultsFrequencies'.
  averageProportion <- c()
  variance <- c()
  
  for (test in unique(allResultsProportions$dataToAnalyse)) {
    df <- allResultsProportions[allResultsProportions$dataToAnalyse==test,]
    
    for (level in unique(df$Expression)) {
      df1 <- df[df$Expression == level,]
      
      if (nrow(df1) >= 1) {
        for (mod in unique(df1$Mod.TF)) {
          df2 <- df1[df1$Mod.TF==mod,]
          
          for (r in unique(df2$Region)) {
            df3 <- df2[df2$Region==r,]
            
            averageProportion <- append(averageProportion, mean(df3$Proportion))
            variance <- append(variance, paste("±", signif(sd(df3$Proportion), digits = 3)))
          }
        }
      } else if (nrow(df1) == 1) {
        averageProportion <- averageProportion
        variance <- variance 
      }
    }
  }
  allResultsFrequencies <- cbind(allResultsFrequencies, data.frame(Mean.Enrichment = averageProportion,
                                                                   Enrichment.Variance = variance))

  if (normalised == FALSE) {
    write.csv(allResultsFrequencies, paste(analysis, "\\Non-normalised\\allResultsFrequencies.csv", sep = "")) 
    write.csv(allResultsProportions, paste(analysis, "\\Non-normalised\\allResultsProportions.csv", sep = "")) 
    
  } else if (normalised == TRUE) {
    write.csv(allResultsFrequencies, paste(analysis, "\\Normalised\\allResultsFrequencies.csv", sep = "")) 
    write.csv(allResultsProportions, paste(analysis, "\\Normalised\\allResultsProportions.csv", sep = "")) 
  }
}

analysisComplete <- TRUE