source("Functions/getGeneCoordinates.R")
source("Functions/Gene width&range.R")
source("Functions/AxisGroup column.R")

# Import genome data.
genomeData <- as.data.frame(read.csv("Wheat data/Wheat protein coding genes.csv"))
genomeData <- genomeData[,-1]

# Import ChIP-seq data.
source("Wheat NLRs/getChIP_seq_data.R")
nextflowOutput <- getChIP_seq_data()

# Convert 'nextflowOutput'dataset into bedfile format.
nextflowOutputBed <- GRanges(nextflowOutput[,c(1:3,8)])

# Import sample gene sets. Store in a list.
sampleGenes <- list()

# Import R-genes and get gene coordinates.
sampleGenes[["R-genes"]] <- as.data.frame(read.xlsx("Wheat NLRs/Wheat NLRs (Andersen et al., 2020).xlsx"))
sampleGenes[["R-genes"]]$Gene <- str_match(sampleGenes[["R-genes"]]$Gene, "^([a-zA-Z0-9]+.*).[0-9]+$")[,2]
  
sampleGenes[["R-genes"]] <- genomeData[which(genomeData$Gene %in% sampleGenes[["R-genes"]]$Gene),]

bedFile <- GRanges(Gene = sampleGenes[["R-genes"]]$Gene,
                   seqnames = sampleGenes[["R-genes"]]$seqnames,
                   ranges = sampleGenes[["R-genes"]]$ranges)

export.bed(bedFile, "Wheat NLRs/Wheat_R-genes.bed")
rm(bedFile)

# Import control gene sets. Store in a list.
# Sample 2000 random control genes normalised for gene length and expression level.
controlGenes <- c()

for (n in seq(from = .1, to = 1, by = .1)) {
  controlGenes <- append(controlGenes, sample(which(genomeData$width > quantile(sampleGenes[["R-genes"]]$width, probs = n-.1) &
                                                      genomeData$width <= quantile(sampleGenes[["R-genes"]]$width, probs = n) &
                                                      genomeData$TPM < max(sampleGenes[["R-genes"]]$TPM)),
                                              nrow(sampleGenes[["R-genes"]])*.1))
}

# Import control gene sets, n Store in a list.
sampleGenes[["Controls"]] <- genomeData[controlGenes,]

# Add 'expressionLevel' column to each dataset.
for (set in names(sampleGenes)) {
  sampleGenes[[set]]$expressionLevel <- rep(NA, times = nrow(sampleGenes[[set]]))
  
  for (row in 1:nrow(sampleGenes[[set]])) {
   if (sampleGenes[[set]][row, "TPM"] <= 1) {
        sampleGenes[[set]][row, "expressionLevel"] <- "No Exp."
   }
    else if (sampleGenes[[set]][row, "TPM"] > 1) {
        sampleGenes[[set]][row, "expressionLevel"] <- "Low Exp."
    }
  }
}

rm(row, set, n)

# Perform enrichment analysis on each gene set.
for (geneSet in names(sampleGenes)) {
  for (level in c("No Exp.", "Low Exp.")) {
    print(paste(geneSet, level))
    
    # Create workbook into which the output data will be saved.
    allProportions_wb <- createWorkbook()
    allFrequency_wb <- createWorkbook()
    
    df <- sampleGenes[[geneSet]][which(sampleGenes[[geneSet]]$expressionLevel==level),]
    
    
    # Get the coordinates for each genomic region.
    geneRegions <- getGeneCoordinates(df, genomeData)
    
    # Use bedtools intersect to find the overlap between ChIP-seq peaks and each gene region.
    for (mod in unique(nextflowOutput$Mod.TF)) {
      addWorksheet(allProportions_wb, sheetName = mod)
      addWorksheet(allFrequency_wb, sheetName = mod)
      
      allProportions_data <- data.frame()
      allFrequency_data <- data.frame()
      
      peaksPerModification <- nextflowOutputBed[which(nextflowOutputBed$Mod.TF==mod),]
      
      for (region in names(geneRegions)) {
        queryBed <- GRanges(geneRegions[[region]][,c("Gene","seqnames","start","end","width")])
        
        overlap <- bt.intersect(peaksPerModification, queryBed, wo = TRUE) 
        
        # Merge multiple peaks overlapping the same region.
        # Determine the proportion of mergedOverlap, with the maximum being 100%.
        
        mergedOverlap <- data.frame()
        
        for (gene in unique(geneRegions[[region]]$Gene)) {
          row <- which(geneRegions[[region]]$Gene == gene)
          
          if (gene %in% overlap$V12) {
            peakOverlaps <- overlap[overlap$V12==gene,]
            
            mergedOverlap <- rbind(mergedOverlap, data.frame(Gene = gene,
                                                             seqnames = geneRegions[[region]][row, "seqnames"],
                                                             start = geneRegions[[region]][row, "start"],
                                                             end = geneRegions[[region]][row, "end"],
                                                             width = geneRegions[[region]][row, "width"],
                                                             overlap = sum(peakOverlaps$V13)/peakOverlaps$V10[1]))
            
          } else mergedOverlap <- rbind(mergedOverlap, data.frame(Gene = gene,
                                                                  seqnames = geneRegions[[region]][row, "seqnames"],
                                                                  start = geneRegions[[region]][row, "start"],
                                                                  end = geneRegions[[region]][row, "end"],
                                                                  width = geneRegions[[region]][row, "width"],
                                                                  overlap = 0))
        }
        
        mergedOverlap[which(mergedOverlap$overlap > 1), "overlap"] <- 1
        
        allProportions_data <- rbind(allProportions_data, data.frame(Gene = mergedOverlap$Gene,
                                                                     seqnames = mergedOverlap$seqnames,
                                                                     start = mergedOverlap$start,
                                                                     end = mergedOverlap$end, 
                                                                     width = mergedOverlap$width,
                                                                     region = rep(region, times = nrow(mergedOverlap)),
                                                                     Mod.TF = rep(mod, times = nrow(mergedOverlap)),
                                                                     overlap = mergedOverlap$overlap,
                                                                     geneSet = rep(paste(geneSet, level), times = nrow(mergedOverlap)),
                                                                     axisGroup = rep(geneRegionAxisLocations(region), times = nrow(mergedOverlap))))
        
        
        allFrequency_data <- rbind(allFrequency_data, data.frame(region = region,
                                                                 Mod.TF = mod,
                                                                 mean.overlap = mean(mergedOverlap$overlap),
                                                                 frequency = signif((nrow(mergedOverlap[which(mergedOverlap$overlap>0.3),])/nrow(mergedOverlap))*100, digits = 3),
                                                                 geneSet = paste(geneSet, level),
                                                                 axisGroup = geneRegionAxisLocations(region)))
        
      }
      writeData(allProportions_wb, sheet = mod, allProportions_data)
      writeData(allFrequency_wb, sheet = mod, allFrequency_data)
    }
    # Save workbook as an .xlsx file.
    saveWorkbook(allProportions_wb, paste("Wheat NLRs/Data/Enrichment results/All proportions/", geneSet, level, ".xlsx", sep = ""), overwrite = TRUE) 
    saveWorkbook(allFrequency_wb, paste("Wheat NLRs/Data/Enrichment results/All frequencies/", geneSet, level, ".xlsx", sep = ""), overwrite = TRUE)  
  }
}

rm(df, nextflowOutputBed, queryBed, allProportions_wb, allProportions_data, allFrequency_wb, allFrequency_data, overlap, region, geneSet, level)

# Plot results.
axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "TTS", "Downstream \n(200bp)", "Intergenic")

# Import results data, merging all gene sets for each genotype. Plot results.
for (mod in unique(nextflowOutput$Mod.TF)) {
  
  allProportions_results <- data.frame()
  allFrequency_results <- data.frame()
  
  for (file in list.files("Wheat NLRs/Data/Enrichment results/All proportions")) {
    allProportions <- read.xlsx(paste("Wheat NLRs/Data/Enrichment results/All proportions/", file, sep = ""), sheet = mod)
    allProportions_results <- rbind(allProportions_results, allProportions)
  }
  for (file in list.files("Wheat NLRs/Data/Enrichment results/All frequencies")) {
    allFrequency <- read.xlsx(paste("Wheat NLRs/Data/Enrichment results/All frequencies/", file, sep = ""), sheet = mod)
    allFrequency_results <- rbind(allFrequency_results, allFrequency)
  }
      
  my_comparisons <- list(c("Controls No Exp.", "R-genes No Exp."), 
                         c("Controls Low Exp.", "R-genes Low Exp."),
                         c("R-genes No Exp.", "R-genes Low Exp."))
  
  stat.test <- allProportions_results %>% group_by(axisGroup) %>% 
    t_test(overlap ~ geneSet, comparisons = my_comparisons) %>% 
    mutate(y.position = rep(c(0.92, 0.92, 0.97), times = 10))
  
  allFrequency_results$geneSet <- factor(allFrequency_results$geneSet, levels = c("Controls No Exp.", "R-genes No Exp.",
                                                                                  "Controls Low Exp.", "R-genes Low Exp."))
  
  plot <- ggbarplot(allFrequency_results, x = "geneSet", y="mean.overlap", ylab = "Average Enrichment",
                    color = "black", fill = "geneSet", lab.size = 12,
                    palette = c("azure3", "azure4", "darkslategray4", "darkslategrey"), 
                    title = mod) + theme_bw() +
    stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", size = 4,
      tip.length = 0.01, hide.ns = FALSE) +
    coord_cartesian(ylim= c(0,1), clip = "off") +
    
    geom_text(data = allFrequency_results, aes(x = geneSet, y = mean.overlap+0.025, label = frequency), size = 2) +
    
    font("title", size = 14) +
    font("ylab", size = 12) +
    font("legend.title", size = 12) +
    font("legend.text", size = 10) +
    font("caption", size = 12) +
    font("axis.text", size = 14)
  
  plot <- facet(plot, facet.by = "axisGroup", nrow = 1, panel.labs.font = list(size = 10),
                panel.labs = list(axisGroup = c("Intergenic","Promotor \n(1kb)","Promotor \n(500bp)", "20%",            
                                                "40%","60%","80%","100%","Downstream \n(200bp)","Intergenic")))
  
  
  plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                font.ytickslab = 8) + theme(legend.key.size = unit(.5, 'cm')) + guides(fill=guide_legend(nrow=1,byrow=TRUE))
  
  pdf(paste("Wheat NLRs/Graphs/Enrichment per region/", mod, ".pdf", sep = ""), width = 12, height = 5)
  print(plot)
  dev.off()
}
