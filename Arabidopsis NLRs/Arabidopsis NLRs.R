source("Functions/getGeneCoordinates.R")
source("Functions/Gene width&range.R")
source("Functions/AxisGroup column.R")

# Import ChIP-seq data.
source("Arabidopsis NLRs/getChIP_seq_data.R")
nextflowOutput <- getChIP_seq_data()

# Convert 'nextflowOutput dataset into bedfile format.
nextflowOutputBed <- GRanges(nextflowOutput[,c(1:3,8)])

# Import genome data.
genomeData <- as.data.frame(read.csv("Arabidopsis NLRs/PlantExp data/geneExpression.csv"))
genomeData <- genomeData[,-1]

# Import sample gene sets. Store in a list.
NLRs <- as.data.frame(read.xlsx("Arabidopsis NLRs/Data/Arabidopsis NLRs.xlsx"))
sampleGenes <- list(NLRs = genomeData[which(genomeData$Gene %in% NLRs$Gene),])

# Import control gene sets. Store in a list.
# Sample ~1000 random genes with an equal distribution of gene lengths, expression levels, 
# and euchromatic/heterochromatic genes as the NLRs.
euchromaticGenes <- c()
heterochromaticGenes <- c()

for (n in seq(from = .1, to = 1, by = .1)) {
  euchromaticGenes <- append(euchromaticGenes, sample(genomeData[which(genomeData$width > quantile(sampleGenes[["NLRs"]]$width, probs = n-.1) &
                                                                           genomeData$width <= quantile(sampleGenes[["NLRs"]]$width, probs = n) &
                                                                           genomeData$Chromosomal.Region=="euchromatic" &
                                                                           !(genomeData$Gene %in% sampleGenes[["NLRs"]]$Gene) &
                                                                           genomeData$TPM < max(sampleGenes[["NLRs"]]$TPM)),"Gene"], 
                                                      100*(nrow(sampleGenes[["NLRs"]][which(sampleGenes[["NLRs"]]$width > quantile(sampleGenes[["NLRs"]]$width, probs = n-.1) &
                                                                                 sampleGenes[["NLRs"]]$width <= quantile(sampleGenes[["NLRs"]]$width, probs = n) &
                                                                                 sampleGenes[["NLRs"]]$Chromosomal.Region=="euchromatic"),])/(nrow(sampleGenes[["NLRs"]])*.1))))
  
  heterochromaticGenes <- append(heterochromaticGenes, sample(genomeData[which(genomeData$width > quantile(sampleGenes[["NLRs"]]$width, probs = n-.1) &
                                                                                 genomeData$width <= quantile(sampleGenes[["NLRs"]]$width, probs = n) &
                                                                                 genomeData$Chromosomal.Region=="heterochromatic" &
                                                                                 !(genomeData$Gene %in% sampleGenes[["NLRs"]]$Gene) &
                                                                                 genomeData$TPM < max(sampleGenes[["NLRs"]]$TPM)),"Gene"], 
                                                              100*(nrow(sampleGenes[["NLRs"]][which(sampleGenes[["NLRs"]]$width > quantile(sampleGenes[["NLRs"]]$width, probs = n-.1) &
                                                                                                      sampleGenes[["NLRs"]]$width <= quantile(sampleGenes[["NLRs"]]$width, probs = n) &
                                                                                                      sampleGenes[["NLRs"]]$Chromosomal.Region=="heterochromatic"),])/(nrow(sampleGenes[["NLRs"]])*.1))))
}
sampleGenes[["Controls"]] <- genomeData[which(genomeData$Gene %in% c(euchromaticGenes, heterochromaticGenes)),]
rm(euchromaticGenes, heterochromaticGenes, NLRs, n)

# Add 'expressionLevel' column to each dataset.
for (set in names(sampleGenes)) {
  sampleGenes[[set]]$expressionLevel <- rep(NA, times = nrow(sampleGenes[[set]]))
  
  for (row in 1:nrow(sampleGenes[[set]])) {
    if (sampleGenes[[set]][row, "TPM"] <= 5) {
      sampleGenes[[set]][row, "expressionLevel"] <- "No Expression"
    }
    else if (sampleGenes[[set]][row, "TPM"] > 5) {
      sampleGenes[[set]][row, "expressionLevel"] <- "Low Expression"
    }
  }
}

rm(row, set)

# Perform enrichment analysis on each gene set.
for (geneSet in names(sampleGenes)) {
  for (level in c("No Expression", "Low Expression")) {
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
    saveWorkbook(allProportions_wb, paste("Arabidopsis NLRs/Data/Enrichment results/All proportions/", geneSet, "_", level, ".xlsx", sep = ""), overwrite = TRUE) 
    saveWorkbook(allFrequency_wb, paste("Arabidopsis NLRs/Data/Enrichment results/All frequencies/", geneSet, "_", level, ".xlsx", sep = ""), overwrite = TRUE) 
  }
}
rm(df, nextflowOutputBed, queryBed, allProportions_wb, allProportions_data, allFrequency_wb, allFrequency_data, overlap, region, geneSet)

# Plot results.
axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "TTS", "Downstream \n(200bp)", "Intergenic")

# Import results data, merging all gene sets for each set. Plot results.
for (mod in unique(nextflowOutput$Mod.TF)) {
  
  allProportions_results <- data.frame()
  allFrequency_results <- data.frame()
  
  for (file in list.files("Arabidopsis NLRs/Data/Enrichment results/All proportions/")) {
    
    allProportions <- read.xlsx(paste("Arabidopsis NLRs/Data/Enrichment results/All proportions/", file, sep = ""), sheet = mod)
    allProportions_results <- rbind(allProportions_results, allProportions)
  }
  for (file in list.files("Arabidopsis NLRs/Data/Enrichment results/All frequencies/")) {
    allFrequency <- read.xlsx(paste("Arabidopsis NLRs/Data/Enrichment results/All frequencies/", file, sep = ""), sheet = mod)
    allFrequency_results <- rbind(allFrequency_results, allFrequency)
  }
  
  my_comparisons <- list(c("Controls No Expression", "NLRs No Expression"), 
                         c("Controls Low Expression", "NLRs Low Expression"),
                         c("NLRs No Expression", "NLRs Low Expression"))
  
  stat.test <- allProportions_results %>% group_by(axisGroup) %>% 
    t_test(overlap ~ geneSet, comparisons = my_comparisons) %>% 
    mutate(y.position = rep(c(0.98, 0.98, 1.06), times = 10))
  
  allFrequency_results$geneSet <- factor(allFrequency_results$geneSet, levels = c("Controls No Expression", "NLRs No Expression",
                                                                                  "Controls Low Expression", "NLRs Low Expression"))
  
  plot <- ggbarplot(allFrequency_results, x = "geneSet", y="mean.overlap", ylab = "Average enrichment",
                    color = "black", fill = "geneSet", 
                    palette = c("azure3", "cadetblue", "bisque2", "lightsalmon2"), 
                    title = mod) + theme_bw() +
    stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", size = 4,
      tip.length = 0.01, hide.ns = FALSE) +
    coord_cartesian(ylim= c(0,1.08), clip = "off") +
    
    geom_text(data = allFrequency_results, aes(x = geneSet, y = mean.overlap+0.025, label = frequency), size = 2) +
    
    font("title", size = 16) +
    font("ylab", size = 14) +
    font("legend.title", size = 12) +
    font("legend.text", size = 10) +
    font("caption", size = 12) 
  
  plot <- facet(plot, facet.by = "axisGroup", nrow = 1, panel.labs.font = list(size = 10),
                panel.labs = list(axisGroup = c("Intergenic","Promotor \n(1kb)","Promotor \n(500bp)", "20%",            
                                                "40%","60%","80%","100%","Downstream \n(200bp)","Intergenic")))
  
  
  plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                font.ytickslab = 8)
  
  pdf(paste("Arabidopsis NLRs/Graphs/Enrichment per region/", mod, ".pdf", sep = ""), width = 12, height = 5)
  print(plot)
  dev.off() 
}
