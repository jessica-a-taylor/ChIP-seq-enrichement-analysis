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
library(ggpubr)
library(rstatix)
library(tidyverse)
library(glue)


# Enrichment analysis based on the occurrence of significant peaks.
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")
source("Functions\\Regions column.R")

ChIP_experiments <- as.data.frame(read_csv("Nextflow_backup/ChIP experiment SRA data.csv"))

# Combine nextflow pipeline outputs into a single dataframe.
nextflowOutput <- data.frame()

for (file in list.files(path = "Nextflow_backup", pattern = "Peaks.bed")) {
  # Merge all broad and narrow peaks datasets.
  data <-  as.data.frame(import.bed(paste("Nextflow_backup/", file, sep = "")))
  
  # Add a column with the experiment code.
  data$experiment <- rep(str_match(file, "^(SRR[0-9]+).*$")[,-1], times = nrow(data))
  
  # Add data to 'nextflowOutput'.
  nextflowOutput <- rbind(nextflowOutput, data)
}

# Change format of 'seqnames' column.
nextflowOutput$seqnames <- str_match(nextflowOutput[,"seqnames"], "^Chr(.)$")[,-1]

# Add ranges column to 'nextflowOutput'.
nextflowOutput$ranges <- mergeCoordinates(nextflowOutput) 

# Add column containing the chromatin modification/TF investigated.
nextflowOutput$`Mod.TF` <- rep(NA, times = nrow(nextflowOutput))

for (mod in unique(ChIP_experiments$`Modification/TF`)) {
  focusModification <- ChIP_experiments[ChIP_experiments$`Modification/TF`==mod,]
  
  # Create a list into which the ChIP-seq experiments for the focus modification/TF will be stored.
  focusExperiments <- c()
  
  for (row in 1:nrow(focusModification)) {
    focusExperiments <- append(focusExperiments, focusModification[row, "Sample data"])
  }
  # Convert space-separated experiment list into a comma-separated list.
  focusExperiments <- glue_collapse(focusExperiments, ",")
  focusExperiments <- gsub(" ", ",", focusExperiments)
  focusExperiments <- c(strsplit(focusExperiments, ",")[[1]])
  
  # Replace 'NAs' in 'Mod.TF' column with the focus modification/TF
  nextflowOutput[which(nextflowOutput$experiment %in% focusExperiments),"Mod.TF"] <- mod
}

rm(focusModification, data, file, mod, n, row)

# Perform enrichment analysis.

# Import coordinates of the genomic regions of interest.
genomicData <- as.data.frame(read_csv("Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Import list of R-genes.
NLR_genes <- as.data.frame(read_xlsx("Arabidopsis NLRs.xlsx", sheet = 1))
NLR_genes <- genomicData[which(genomicData$Gene %in% NLR_genes$Gene),]

# Sample 1000 random control genes, then sort R-genes and control genes based on expression level.
# Ensure that the number of R-genes and control genes is the same for a particular expression level.
source("Functions\\PlantExp.R")

for (normalised in c(TRUE, FALSE)) {
  sampleGenes <- PlantExp(genomicData, NLR_genes)

  # Plot average gene size between R-genes and control genes.
  geneWidth <- data.frame()
  
  for (set in names(sampleGenes)) {
    geneWidth <- rbind(geneWidth, data.frame(Gene = sampleGenes[[set]]$Gene,
                                             GeneSet = sampleGenes[[set]]$GeneSet,
                                             GeneWidth = sampleGenes[[set]]$width/1000))
    
  }
  
  plot <- ggplot(geneWidth, aes(x = GeneSet, y = GeneWidth)) +
    geom_boxplot() + labs(x = "Gene set", y = "Gene width (kb)") +
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = "R-gene") + theme_bw() 
  
  if (normalised == TRUE) {
    write.csv(geneWidth, paste("PlantExp data\\Normalised\\geneWidth.csv", sep = "")) 
    ggsave("Graphs\\Gene width comparison.png", plot = plot, width = 8, height = 4)  
  } else if (normalised == FALSE) {
    write.csv(geneWidth, paste("PlantExp data\\Non-normalised\\geneWidth.csv", sep = "")) 
    ggsave("Graphs\\Non-normalised\\Gene width comparison.png", plot = plot, width = 4, height = 3)  
  }
  jobRunScript("Script for analysis.R", importEnv = TRUE)
}

for (normalised in c(TRUE, FALSE)) {
  jobRunScript("Plot enrichment.R",  importEnv = TRUE)
  jobRunScript("Compare average gene size.R",  importEnv = TRUE) 
}