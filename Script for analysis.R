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


if (analysis == "PlantExp data") {
  dataToAnalyse <- sampleGenesPlantExp
} else if (analysis == "RNA-seq data") {
  dataToAnalyse <- sampleGenesRNAseq
}

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
    
    # Determine the proportion of each gene region overlapping with a significant peak.
    geneRegions <- getGeneCoordinates(geneSet)
    
    proportionPerRegion <- proportionsFunction(geneRegions, allOverlaps, nextflowOutput)
    
    # Add a column to 'proportionPerRegion' with the numbers for 
    # each gene region that will correspond with their position on the x axis.
    proportionPerRegion <- geneRegionAxisLocations(proportionPerRegion, geneRegions)
    
    rm(geneRegions)
    
    # Add a column to 'proportionPerRegion' with the current expression level.
    proportionPerRegion <- expressionColumn(proportionPerRegion, test)
    
    # Store final results in 'dataToAnalyseProportions'.
    dataToAnalyseProportions[[test]] <- proportionPerRegion
    
    print(test)
  }
  
  # Merge all data from all sample gene sets into one big dataframe.
  allResultsProportions <- data.frame()
  
  for (test in names(dataToAnalyseProportions)) {
    
    df <- dataToAnalyseProportions[[test]]
    df <- cbind(df, data.frame(dataToAnalyse = rep(unlist(str_split(test, "_"))[1], times = nrow(df))))
    
    allResultsProportions <- rbind(allResultsProportions, df)
  }
  
  # Calculate the mean proportion of overlap and add as a new column to the dataframe.
  allResultsAverageProportions <- data.frame()
  
  for (test in unique(allResultsProportions$dataToAnalyse)) {
    df <- allResultsProportions[allResultsProportions$Sample==test,]
    
    for (level in unique(df$Expression)) {
      df1 <- df[df$Expression == level,]
      
      if (nrow(df1) >= 1) {
        for (exp in unique(df1$Experiment)) {
          df2 <- df1[df1$Experiment==exp,]
          
          for (r in unique(df2$Region)) {
            df3 <- df2[df2$Region==r,]
            
            allResultsAverageProportions <- rbind(allResultsAverageProportions, data.frame(Region = r,
                                                                                           Experiment = exp,
                                                                                           Proportion = mean(df3$Proportion),
                                                                                           axisGroup = df3$axisGroup[1],
                                                                                           Expression = level,
                                                                                           Sample = test,
                                                                                           SampleSize = nrow(df3)))
          }
        }
      } else allResultsAverageProportions <- allResultsAverageProportions
    }
  }
  
  
  # Average enrichment across all genes at all expression levels.
  proportionsAllExpressionLevels <- data.frame()
  
  for (test in unique(allResultsProportions$dataToAnalyse)) {
    df1 <- allResultsProportions[allResultsProportions$dataToAnalyse==test,]
    
    if (nrow(df1) >= 1) {
      for (exp in unique(df1$Experiment)) {
        df2 <- df1[df1$Experiment==exp,]
        
        for (r in unique(df2$Region)) {
          df3 <- df2[df2$Region==r,]
          
          proportionsAllExpressionLevels <- rbind(proportionsAllExpressionLevels, data.frame(Region = r,
                                                                                             Experiment = exp,
                                                                                             Proportion = mean(df3$Proportion),
                                                                                             axisGroup = df3$axisGroup[1],
                                                                                             Sample = test,
                                                                                             SampleSize = nrow(df3)))
          
        }
      }
    } else proportionsAllExpressionLevels <- proportionsAllExpressionLevels
  }
  if (normalised == FALSE) {
    write.csv(allResultsProportions, paste(analysis, "\\Non-normalised\\allResultsProportions.csv", sep = "")) 
  } else if (normalised == TRUE) {
    write.csv(allResultsProportions, paste(analysis, "\\Normalised\\allResultsProportions.csv", sep = "")) 
  }
}