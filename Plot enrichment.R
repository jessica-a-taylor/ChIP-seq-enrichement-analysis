library(ggplot2)
library(data.table)
library(grid)
library(readr)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(rstudioapi)

axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "TTS", "Downstream \n(200bp)", "Intergenic")

if (normalised == TRUE) {
  allFrequencies <- data.frame(read_csv(paste("PlantExp data\\Normalised\\allFrequencies.csv", sep = "")))
  promoterEnrichment <- data.frame(read_csv(paste("PlantExp data\\Normalised\\promoterEnrichment.csv", sep = "")))  
  genebodyEnrichment <- data.frame(read_csv(paste("PlantExp data\\Normalised\\genebodyEnrichment.csv", sep = "")))
  allProportionsPerRegion <- data.frame(read_csv(paste("PlantExp data\\Normalised\\allProportionsPerRegion.csv", sep = "")))  
  allProportionsPerRegion <- allProportionsPerRegion[,-1]
  
} else if (normalised == FALSE) {
  allFrequencies <- data.frame(read_csv(paste("PlantExp data\\Non-normalised\\allFrequencies.csv", sep = "")))
  promoterEnrichment <- data.frame(read_csv(paste("PlantExp data\\Non-normalised\\promoterEnrichment.csv", sep = "")))  
  genebodyEnrichment <- data.frame(read_csv(paste("PlantExp data\\Non-normalised\\genebodyEnrichment.csv", sep = "")))
  allProportionsPerRegion <- data.frame(read_csv(paste("PlantExp data\\Non-normalised\\allProportionsPerRegion.csv", sep = "")))  
  allProportionsPerRegion <- allProportionsPerRegion[,-1]
}
  
# Replace comma in 'Comparisons' column with \n.
allFrequencies$GeneSet <- paste(str_match(allFrequencies$GeneSet, "^([A-Za-z]+.gene).*$")[,-1], " \n", 
                                       str_match(allFrequencies$GeneSet, "^[A-Za-z]+.gene(.*)$")[,-1], sep = "")

allProportionsPerRegion$GeneSet <- paste(str_match(allProportionsPerRegion$GeneSet, "^([A-Za-z]+.gene).*$")[,-1], " \n", 
                                       str_match(allProportionsPerRegion$GeneSet, "^[A-Za-z]+.gene(.*)$")[,-1], sep = "")

allFrequencies <- allFrequencies[order(factor(allFrequencies$GeneSet, 
                                                            levels = c("Control gene \n No Expression", "R-gene \n No Expression",
                                                                       "Control gene \n Low Expression", "R-gene \n Low Expression"))),]

allProportionsPerRegion <- allProportionsPerRegion[order(factor(allProportionsPerRegion$GeneSet, 
                                                            levels = c("Control gene \n No Expression", "R-gene \n No Expression",
                                                                       "Control gene \n Low Expression", "R-gene \n Low Expression"))),]


my_comparisons <- list(c("Control gene \n No Expression", "R-gene \n No Expression"), 
                       c("Control gene \n Low Expression", "R-gene \n Low Expression"),
                       c("R-gene \n No Expression", "R-gene \n Low Expression"))

# Plot bar graph of average enrichment per region.
for (mod in unique(allProportionsPerRegion$Mod.TF)) {
  df <- allFrequencies[allFrequencies$Mod.TF==mod,]
  
  comparison_df <- allProportionsPerRegion[allProportionsPerRegion$Mod.TF==mod,]
  
  stat.test <- comparison_df %>% group_by(axisGroup) %>% 
    t_test(Proportion ~ GeneSet, comparisons = my_comparisons) %>% 
    mutate(y.position = rep(c(0.98, 0.98, 1.06), times = 10))
  
  plot <- ggbarplot(df, x = "GeneSet", y="Enrichment.mean", ylab = "Average enrichment",
                    color = "black", fill = "GeneSet", 
                    palette = c("azure3", "cadetblue", "bisque2", "lightsalmon2"), 
                    title = mod) + theme_bw() +
    stat_pvalue_manual(
      stat.test, 
      label = "p.adj.signif", size = 4,
      tip.length = 0.01, hide.ns = FALSE) +
    coord_cartesian(ylim= c(0,1.08), clip = "off") +
    
    geom_text(data = df, aes(x = GeneSet, y = Enrichment.mean+0.025, label = Count), size = 2) +
    
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
  
  ggsave(paste("Graphs\\Normalised\\Enrichment per region\\", mod, ".png", sep = ""), plot = plot, width = 10, height = 4)  
}

# Plot bar graph of enrichment per R-gene promoter.
for (level in unique(promoterEnrichment$ExpressionLevel)) {
  df <- promoterEnrichment[promoterEnrichment$ExpressionLevel==level,]
  df <- df[which(df$Mod.TF %in% c("H3K4me3", "H3K36me3", "H3K9ac", "H3K27ac", "H3K27me3", "H2A.Z", "H2Bub", "H2AK121ub", "H3K4me1", "H3K9me2")),]
  
  plot <- ggplot(df, aes(Mod.TF, Gene)) +
    geom_tile(aes(fill = Proportion), colour = "white") +
    scale_fill_gradient(low = "white", high = "indianred2") +
    labs(x = "Histone modification", y = "Gene", fill = "Enrichment") +
    theme(axis.text.y = element_text(size = 11, colour = "grey13"),
          axis.text.x = element_text(size = 13, colour = "grey13"),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 

  ggsave(paste("Graphs\\Normalised\\Enrichment per gene\\",level,"promotor epigenetic enrichment.png", sep = ""), plot = plot, width = 16, height = 16)
  
  df <- promoterEnrichment[promoterEnrichment$ExpressionLevel==level,]
  df <- df[which(df$Mod.TF %in% c("WRKY18", "WRKY33", "WRKY40", "VAL1", "VAL2", "TPR1", "ATXR7", "HDA9", "REF6", "EDM2",
                                  "LPH1", "SPTL6", "GBPL3", "CCA1", "TOC1", "LHY", "NUP1", "CRWN1")),]
  
  plot <- ggplot(df, aes(Mod.TF, Gene)) +
    geom_tile(aes(fill = Proportion), colour = "white") +
    scale_fill_gradient(low = "white", high = "indianred2") +
    labs(x = "Histone modification", y = "Gene", fill = "Enrichment") +
    theme(axis.text.y = element_text(size = 11, colour = "grey13"),
          axis.text.x = element_text(size = 13, colour = "grey13"),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 

  ggsave(paste("Graphs\\Normalised\\Enrichment per gene\\",level, "promotor TF enrichment.png", sep = ""), plot = plot, width = 16, height = 16)
}

# Plot bar graph of enrichment per R-gene body
for (level in unique(genebodyEnrichment$ExpressionLevel)) {
  df <- genebodyEnrichment[genebodyEnrichment$ExpressionLevel==level,]
  df <- df[which(df$Mod.TF %in% c("H3K4me3", "H3K36me3", "H3K9ac", "H3K27ac", "H3K27me3", "H2A.Z", "H2Bub", "H2AK121ub", "H3K4me1", "H3K9me2")),]
  
  plot <- ggplot(df, aes(Mod.TF, Gene)) +
    geom_tile(aes(fill = Proportion), colour = "white") +
    scale_fill_gradient(low = "white", high = "indianred2") +
    labs(x = "Histone modification", y = "Gene", fill = "Enrichment") +
    theme(axis.text.y = element_text(size = 11, colour = "grey13"),
          axis.text.x = element_text(size = 13, colour = "grey13"),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 
  
  ggsave(paste("Graphs\\Normalised\\Enrichment per gene\\",level,"genebody epigenetic enrichment.png", sep = ""), plot = plot, width = 16, height = 16)
  
  df <- genebodyEnrichment[genebodyEnrichment$ExpressionLevel==level,]
  df <- df[which(df$Mod.TF %in% c("WRKY18", "WRKY33", "WRKY40", "VAL1", "VAL2", "TPR1", "ATXR7", "HDA9", "REF6", "EDM2",
                                  "LPH1", "SPTL6", "GBPL3", "CCA1", "TOC1", "LHY", "NUP1", "CRWN1")),]
  
  plot <- ggplot(df, aes(Mod.TF, Gene)) +
    geom_tile(aes(fill = Proportion), colour = "white") +
    scale_fill_gradient(low = "white", high = "indianred2") +
    labs(x = "Histone modification", y = "Gene", fill = "Enrichment") +
    theme(axis.text.y = element_text(size = 11, colour = "grey13"),
          axis.text.x = element_text(size = 13, colour = "grey13"),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) 
  
  ggsave(paste("Graphs\\Normalised\\Enrichment per gene\\",level, "genebody TF enrichment.png", sep = ""), plot = plot, width = 16, height = 16)
}