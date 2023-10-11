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

setwd("C:/Users/jexy2/OneDrive/ChIP-seq-enrichment-analysis")

# Enrichment analysis based on the occurrence of significant peaks.
source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")
source("Functions\\Regions column.R")

ChIP_experiments <- as.data.frame(read_xlsx("Nextflow_backup/ChIP experiment SRA data.xlsx"))

# Combine nextflow pipeline outputs into a single dataframe.
nextflowOutput <- data.frame()

for (file in list.files(path = "Nextflow_backup", pattern = "Peaks.bed")) {
  
  # Rename file to include target modification/TF (if it has not been changed already).
  if (str_detect(file, "_") == FALSE) {
    
    file.rename(paste("Nextflow_backup/", file, sep = ""), 
                paste("Nextflow_backup/", ChIP_experiments[which(str_detect(ChIP_experiments$`Sample data`, str_match(file,"^(SRR.*)merged.*$")[,2])==TRUE), "Modification/TF"],
                      "_", file, sep = "")) 
      
    file <- paste(ChIP_experiments[which(str_detect(ChIP_experiments$`Sample data`, str_match(file,"^(SRR.*)merged.*$")[,2])==TRUE), "Modification/TF"],
                    "_", file, sep = "")
  }
  
  # Merge all broad and narrow peaks datasets.
  data <-  as.data.frame(import.bed(paste("Nextflow_backup/", file, sep = "")))
  
  # Add a column with the experiment code.
  data$experiment <- rep(str_match(file, "^.*_(SRR[0-9]+).*$")[,-1], times = nrow(data))
  
  # Add data to 'nextflowOutput'.
  nextflowOutput <- rbind(nextflowOutput, data)
}

# Remove rows for peaks on mitochondrial/chloroplast genome.
nextflowOutput <- nextflowOutput[-which(nextflowOutput$seqnames %in% c("ChrC", "ChrM")),]

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

# Remove rows with any remaining 'NAs'.
nextflowOutput <- nextflowOutput[c(which(!is.na(nextflowOutput$Mod.TF))),]

rm(focusModification, data, file, mod, row)

# Sample 1000 random control genes, then sort R-genes and control genes based on expression level.
# Ensure that the number of R-genes and control genes is the same for a particular expression level.
source("Functions\\PlantExp.R")

sampleGenes <- PlantExp()
  
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

write.csv(geneWidth, paste("PlantExp data\\geneWidth.csv", sep = "")) 
ggsave("Graphs\\Gene width comparison.png", plot = plot, width = 8, height = 4)  


# Perform enrichment analysis.
jobRunScript("Script for analysis.R", importEnv = TRUE)

# Generate enrichment plots.
jobRunScript("Plot enrichment.R",  importEnv = TRUE)
jobRunScript("Compare average gene size.R",  importEnv = TRUE)