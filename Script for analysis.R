library(readxl)
library(karyoploteR)
library(rtracklayer)
library(dplyr)
library(stringr)
library(hash)
library(sets)
library(ggplot2)
library(data.table)
library(grid)
library(readr)
library(rstudioapi)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(glue)
library(TxDb.Athaliana.BioMart.plantsmart28)

# Create workbook into which the output data will be saved.
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

# Save workbook as an .xlsx file.
saveWorkbook(wb, paste(set, ".xlsx", sep = ""), overwrite = TRUE)
