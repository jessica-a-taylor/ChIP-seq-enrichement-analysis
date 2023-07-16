# Load required libraries.

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

# Import coordinates of the genomic regions of interest.
# Options: euchromaticRegions, heterochromaticRegions, euchromaticWithoutTEs

genomicData <- as.data.frame(read_csv("Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Get coordinates for 10 sets of random genes and store in a hash.
source("Functions\\Sample random genes.R")

if (TRUE %in% grepl("control", list.files(path = paste("Data\\RNA-seq data\\"),pattern="*.txt"))) {
  sampleGenes <- existingSets(genomicData)
} else sampleGenes <- geneSets(genomicData)

rm(geneSets, existingSets)

# Import list of R-genes.
NLRgenes <- as.data.frame(read_xlsx("Data\\Arabidopsis NLRs.xlsx", sheet = 1))
NLRgenes <- genomicData[which(genomicData$Gene %in% NLRgenes$Gene),]

# Add R-genes to sampleGenes.
sampleGenes[["NLRs"]] <- NLRgenes

# Get filtered expression data for each set of sample genes. 
# Add dataframes to new sampleGenes hashes for gene sets with particular expression levels.
exLevel <- c("No Expression", "Low Expression", "Intermediate Expression",
             "High Expression")

source("Functions\\PlantExp.R")

sampleGenesPlantExp <- PlantExp(sampleGenes[c("control1","control10","control2","control3","control4",
                                              "control5","control6","control7","control8","control9","NLRs")], 
                                exLevel, "sampleGenes")

source("Functions\\RNA-seq data.R")
sampleGenesRNAseq <- RNA_seqAnalysis(sampleGenes[c("control1","control10","control2","control3","control4",
                                                   "control5","control6","control7","control8","control9","NLRs")], 
                                     exLevel, "sampleGenes")


rm(PlantExp, RNA_seqAnalysis)


source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Coordinates per gene region.R")

for (file in list.files(pattern = "final_counts.csv")) {
  # Import bed file.
  data <- as.data.frame(read_csv(file))
  data <- data[,-c(5:9)]
  
  # Change format of 'seqnames' column.
  data$seqnames <- str_match(data[,"seqnames"], "^Chr(.)$")[,-1]
  
  # Isolate read counts.
  data$counts <- as.numeric(str_match(data[,"counts"], "^[0-9]+.([0-9]+)$")[,-1])
  
  # Add a column with RPK value.
  data$RPK <- data$counts/((data$end - data$start)/1000)
  
  # Calculate sum of RPK values and divide by 1000000 to get the 'per million' scaling factor.
  scalingFactor <- sum(data$RPK)/1000000
  
  # Add a column with normalised read counts (TPM).
  data$TPM <- data$counts/scalingFactor
  
  # Remove fragments from mitochondria and chloroplast genomes.
  data <- data[-c(which(data$seqnames %in% c("M", "C"))),]
  data$seqnames <- as.numeric(data$seqnames)
  
 # Add a column with the range of each fragment.
  data$ranges <- mergeCoordinates(data) 
  
  # Filter for genes of interest.
}