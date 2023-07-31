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
library(TxDb.Athaliana.BioMart.plantsmart28)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(glue)

# Import coordinates of the genomic regions of interest.
# Options: euchromaticRegions, heterochromaticRegions, euchromaticWithoutTEs

genomicData <- as.data.frame(read_csv("Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Get the coordinates for the regions of interest for each gene.
source("Functions\\Coordinates per gene region.R")
geneRegions <- getGeneCoordinates(genomicData)

# Import list of R-genes.
NLR_genes <- as.data.frame(read_xlsx("Arabidopsis NLRs.xlsx", sheet = 1))
NLR_genes <- genomicData[which(genomicData$Gene %in% NLR_genes$Gene),]

# Sample 1000 random control genes, then sort R-genes and control genes based on expression level.
# Ensure that the number of R-genes and control genes is the same for a particular expression level.
source("Functions\\PlantExp.R")
sampleGenesPlantExp <- PlantExp(genomicData, NLR_genes)

rm(PlantExp)


source("Functions\\Get range - merge gene coordinates.R")
source("Functions\\Expression column.R")
source("Functions\\Regions column.R")

ChIP_experiments <- as.data.frame(read_csv("Nextflow_backup/ChIP experiment SRA data.csv"))

axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "TTS", "Downstream \n(200bp)", "Intergenic")


# Enrichment analysis based on the occurrence of significant peaks.
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

# Perform enrichment analysis.
analysis <- "PlantExp data"
jobRunScript("Script for analysis.R", importEnv = TRUE)

jobRunScript("Plot enrichment.R",  importEnv = TRUE)
jobRunScript("Compare average gene size.R", importEnv = TRUE)



# Enrichment analysis based on normalised read counts.
# Perform enrichment analysis.
for (analysis in c("PlantExp data", "RNA-seq data")) {
  if (analysis == "PlantExp data") {
    dataToAnalyse <- sampleGenesPlantExp
  } else if (analysis == "RNA-seq data") {
    dataToAnalyse <- sampleGenesRNAseq
  }
  
  # Create an empty hash into which the results for each experiment will be stored.
  allExperiments <- hash()
  
  for (file in list.files(path = "Nextflow_backup", pattern = "final_counts.csv")) {
    # Import bed file.
    data <- as.data.frame(read_csv(paste("Nextflow_backup/", file, sep = "")))
    data <- data[,-c(5:9)]
    
    # Change format of 'seqnames' column.
    data$seqnames <- str_match(data[,"seqnames"], "^Chr(.)$")[,-1]
    
    # Remove fragments from mitochondria and chloroplast genomes.
    data <- data[-c(which(data$seqnames %in% c("M", "C"))),]
    data$seqnames <- as.numeric(data$seqnames)
    
    # Add a column with the range of each fragment.
    data$ranges <- mergeCoordinates(data) 
    
    # Isolate read counts.
    data$counts <- as.numeric(str_match(data[,"counts"], "^[0-9]+.([0-9]+)$")[,-1])
    
    # Create an empty dataframe into which the total RPK values for each region will be stored.
    totalRPK <- data.frame()
    
    # For each gene in each region...
    for (region in names(geneRegions)) {
      RPK <- c()
      
      for (row in 1:nrow(geneRegions[[region]])) {
        
        # Find peaks overlapping with each region of each gene.
        overlaps <- which(data[,"start"] > as.numeric(geneRegions[[region]][row,"start"]) &
                          data[,"end"] < as.numeric(geneRegions[[region]][row,"end"]) &
                          data[,"seqnames"] == geneRegions[[region]][row,"seqnames"])
        
        # Calculate RPK (total overlapping reads divided by the with of the gene region in kb).
        if (length(overlaps) == 0) {
          RPK <- append(RPK, 0)
        } else if (length(overlaps) >= 1) {
          RPK <- append(RPK, sum(data[c(overlaps),"counts"])/(geneRegions[[region]][row,"width"]/1000))
        }
      }
      totalRPK <- rbind(totalRPK, data.frame(Region = region,
                                             RPK = sum(RPK)))
      
      # Add the RPK, TPM and logTPM values to the 'geneRegions[[region]]' dataframe.
      geneRegions[[region]]$RPK <- RPK
      
      TPM <- c()
      for (row in 1:nrow(geneRegions[[region]])) {
        if (geneRegions[[region]][row,"RPK"] == 0) {
          TPM <- append(TPM, 0)
        } else if (geneRegions[[region]][row,"RPK"] > 0) {
          TPM <- append(TPM, geneRegions[[region]][row,"RPK"]/(sum(totalRPK$RPK)/1000000))
        }
      }
      geneRegions[[region]]$TPM <- TPM
      
      logTPM <- c()
      for (row in 1:nrow(geneRegions[[region]])) {
        if (geneRegions[[region]][row,"TPM"] == 0) {
          logTPM <- append(logTPM, 0)
        } else if (geneRegions[[region]][row,"TPM"] > 0) {
          logTPM <- append(logTPM, log(geneRegions[[region]][row,"TPM"]))
        }
      }
      geneRegions[[region]]$logTPM <- logTPM
    }
    
    # Create an empty dataframe into which the enrichment data for each gene set will be combined.
    allGenes <- data.frame()
    
    for (test in names(dataToAnalyse)[grepl("_", names(dataToAnalyse))]) {
      geneSet <- dataToAnalyse[[test]]
      
      # Create an empty dataframe into which the enrichment data for each region will be combined.
      allRegions <- data.frame()
      
      for (region in names(geneRegions)) {
        # Filter for genes in 'geneSet'.
        genesOfInterest <- geneRegions[[region]][c(which(geneRegions[[region]]$Gene %in% geneSet$Gene)),]
        
        # Add a column with the region of interest.
        genesOfInterest$Region <- rep(region, times = nrow(genesOfInterest))
        
        # Add a column with expression level.
        genesOfInterest <- expressionColumn(genesOfInterest, test)
        
        # Combine data for each gene region into 'allRegions' dataframe.
        allRegions <- rbind(allRegions, data.frame(Gene = genesOfInterest$Gene,
                                                   seqnames = genesOfInterest$seqnames,
                                                   strand = genesOfInterest$strand,
                                                   start = genesOfInterest$start,
                                                   end = genesOfInterest$end,
                                                   width = genesOfInterest$width,
                                                   ranges = genesOfInterest$ranges,
                                                   Region = genesOfInterest$Region,
                                                   Expression = genesOfInterest$Expression,
                                                   RPK = genesOfInterest$RPK,
                                                   TPM = genesOfInterest$TPM,
                                                   logTPM = genesOfInterest$logTPM))
      }
      # Add a column with the x axis coordinates corresponding to each region.
      allRegions <- allRegions[order(factor(allRegions$Region, levels = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
                                                                          "Gene20", "Gene40", "Gene60", "Gene80", "Gene100",
                                                                          "Downstream", "DownstreamIntergenic"))),]
      
      allRegions <- geneRegionAxisLocations(allRegions, unique(allRegions$Region))
      
      # Combine data for each gene region into 'allRegions' dataframe.
      allGenes <- rbind(allGenes, allRegions)
    }
    # Store the 'allGenes' dataframe in the 'allExperiments' hash.
    allExperiments[[str_match(file, "^(SRR[0-9]+).*$")[,-1]]] <- allGenes
  }
  rm(allGenes)
  
  # Combine the results of all experiments for the same modification/TF.
  
  for (mod in unique(ChIP_experiments$`Modification/TF`)) {
    # Create an empty dataframe into which the TPM data for each experiment for a 
    # particular modification/TF will be combined.
    combinedExperiments <- data.frame()
    
    focusModification <- ChIP_experiments[ChIP_experiments$`Modification/TF`==mod,]
    
    # Create a list into which the ChIP-seq experiments for the focus modification/TF will be stored.
    focusExperiments <- c()
    
    for (row in 1:nrow(focusModification)) {
      focusExperiments <- append(focusExperiments, focusModification[row, "Sample data"])
    }
    focusExperiments <- glue_collapse(focusExperiments, ",")
    focusExperiments <- gsub(" ", ",", focusExperiments)
    focusExperiments <- c(strsplit(focusExperiments, ",")[[1]])
    
    # Combine the data from each experiment.
    for (exp in focusExperiments) {
      combinedExperiments <- rbind(combinedExperiments, allExperiments[[exp]])
    }
    
    # Average the TPM value for each gene between experiments.
    combinedAverage <- data.frame()
    
    for (region in unique(combinedExperiments$Region)) {
      df1 <- combinedExperiments[combinedExperiments$Region==region,]
      
      for (gene in unique(df1$Gene)) {
        df2 <- df1[df1$Gene==gene,]
        
        combinedAverage <- rbind(combinedAverage, data.frame(Gene = df2[1, "Gene"],
                                                             Region = df2[1, "Region"],
                                                             Expression = df2[1, "Expression"],
                                                             logTPM = mean(df2[, "logTPM"]),
                                                             axisGroup = df2[1, "axisGroup"]))
      }
    }
    
    # Plot bar graph.
    if (nrow(combinedAverage) >= 1) {
      controlGenes <- combinedAverage[c(which(grepl("control", combinedAverage$Expression) & grepl("No Expression", combinedAverage$Expression))),]
      controlGenes <- rbind(controlGenes, combinedAverage[c(which(grepl("control", combinedAverage$Expression) & grepl("Low Expression", combinedAverage$Expression))),])
      controlGenes$Sample <- rep("Control gene", times = nrow(controlGenes))
      controlGenes$Expression <- str_match(controlGenes$Expression, "^control.*_(.*)$")[,-1]
      
      Rgenes <- combinedAverage[c(which(grepl("NLR", combinedAverage$Expression) & grepl("No Expression", combinedAverage$Expression))),]
      Rgenes <- rbind(Rgenes, combinedAverage[c(which(grepl("NLR", combinedAverage$Expression) & grepl("Low Expression", combinedAverage$Expression))),])
      Rgenes$Sample <- rep("R-gene", times = nrow(Rgenes))
      Rgenes$Expression <- str_match(Rgenes$Expression, "^NLR.*_(.*)$")[,-1]
      
      allGenes <- controlGenes
      allGenes <- rbind(allGenes, Rgenes)
      
      ExpressionGenes <- data.frame()
      
      for (level in unique(allGenes$Expression)) {
        controlGenes1 <- controlGenes[controlGenes$Expression == level,]
        Rgenes1 <- Rgenes[Rgenes$Expression == level,]
        
        for (region in unique(allGenes$Region)) {
          controlGenes2 <- controlGenes1[controlGenes1$Region == region,]
          Rgenes2 <- Rgenes1[Rgenes1$Region == region,]
          
          ExpressionGenes <- rbind(ExpressionGenes, data.frame(Region = rep(region, time = 2),
                                                               logTPM = c(mean(controlGenes2$logTPM), mean(Rgenes2$logTPM)),
                                                               axisGroup = rep(unique(controlGenes2$axisGroup), times = 2),
                                                               Expression = rep(unique(controlGenes2$Expression), times = 2),
                                                               Sample = c("Control gene", "R-gene")))
        }
      }
      
      ExpressionGenes$Comparison <- paste(ExpressionGenes$Sample, ExpressionGenes$Expression, sep = " \n")
      allGenes$Comparison <- paste(allGenes$Sample, allGenes$Expression, sep = " \n")
      
      ExpressionGenes <- ExpressionGenes[order(factor(ExpressionGenes$Comparison, levels = c("Control gene \nNo Expression", "R-gene \nNo Expression", 
                                                                                             "Control gene \nLow Expression", "R-gene \nLow Expression"))),]
      
      allGenes <- allGenes[order(factor(allGenes$Comparison, levels = c("Control gene \nNo Expression", "R-gene \nNo Expression", 
                                                                        "Control gene \nLow Expression", "R-gene \nLow Expression"))),]
      
      
      my_comparisons <- list(c("Control gene \nNo Expression", "R-gene \nNo Expression"), 
                             c("Control gene \nLow Expression", "R-gene \nLow Expression"), 
                             c("R-gene \nNo Expression", "R-gene \nLow Expression"))
      
      stat.test <- allGenes %>% group_by(axisGroup) %>% 
        t_test(logTPM ~ Comparison, comparisons = my_comparisons) %>% 
        mutate(y.position = rep(c(0.7, 0.7, 0.8), times = 10))
      
      plot <- ggbarplot(ExpressionGenes, x = "Comparison", y="logTPM", ylab = "Average TPM per region",
                        color = "black", fill = "Comparison", 
                        palette = c("azure3", "cadetblue", "bisque2", "darksalmon"), title = mod) + 
        stat_pvalue_manual(
          stat.test, 
          label = "p.adj.signif", size = 4,
          tip.length = 0.01, hide.ns = FALSE) +
        coord_cartesian(ylim= c(0,1), clip = "off") +
        theme_bw() +
        
        font("title", size = 16) +
        font("ylab", size = 14) +
        font("legend.title", size = 14) +
        font("legend.text", size = 12) +
        font("caption", size = 12) 
      
      plot <- facet(plot, facet.by = "axisGroup", nrow = 1, panel.labs.font = list(size = 10),
                    panel.labs = list(axisGroup = c("Intergenic","Promotor \n(1kb)","Promotor \n(500bp)", "20%",            
                                                    "40%","60%","80%","100%","Downstream \n(200bp)","Intergenic")))
      
      plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "Gene set",
                    font.ytickslab = 12)
      
      #ggsave(paste("Graphs\\Enrichment\\nextflow\\", analysis, "\\", mod, "_Controls vs R-genes.png", sep = ""), plot = plot, width = 10, height = 4)   
    }
  }
}
  