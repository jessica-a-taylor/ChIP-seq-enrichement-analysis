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

for (normalised in c(TRUE, FALSE)) {
  if (normalised == FALSE) {
    allResultsFrequencies <- data.frame(read_csv(paste(analysis, "\\Non-normalised\\allResultsFrequencies.csv", sep = "")))
    allResultsProportions <- data.frame(read_csv(paste(analysis, "\\Non-normalised\\allResultsProportions.csv", sep = "")))
  } else if (normalised == TRUE) {
    allResultsFrequencies <- data.frame(read_csv(paste(analysis, "\\Normalised\\allResultsFrequencies.csv", sep = "")))
    allResultsProportions <- data.frame(read_csv(paste(analysis, "\\Normalised\\allResultsProportions.csv", sep = "")))
  }
  
  allResultsFrequencies <- allResultsFrequencies[-c(which(allResultsFrequencies$Expression %in% c("Intermediate Expression", "High Expression"))),]
  
  # Plot bar graph.
  ExpressionGenes <- data.frame()
  
  allGenes <- allResultsProportions[which(allResultsProportions$Expression == "No Expression" |
                                            allResultsProportions$Expression == "Low Expression"),]
  
  controlGenes <- allGenes[grepl("control", allGenes$dataToAnalyse),]
  controlGenes$dataToAnalyse <- rep("Control gene", times = nrow(controlGenes))
  
  Rgenes <- allGenes[grepl("NLR", allGenes$dataToAnalyse),]
  Rgenes$dataToAnalyse <- rep("R-gene", times = nrow(Rgenes))
  
  allGenes <- controlGenes
  allGenes <- rbind(allGenes, Rgenes)
  
  for (level in unique(allGenes$Expression)) {
    controlGenes1 <- controlGenes[controlGenes$Expression == level,]
    Rgenes1 <- Rgenes[Rgenes$Expression == level,]
    
    for (mod in unique(allGenes$`Mod.TF`)) {
      
      controlGenes2 <- controlGenes1[controlGenes1$`Mod.TF` == mod,]
      Rgenes2 <- Rgenes1[Rgenes1$`Mod.TF` == mod,]
      
      for (region in unique(allGenes$Region)) {
        controlGenes3 <- controlGenes2[controlGenes2$Region == region,]
        Rgenes3 <- Rgenes2[Rgenes2$Region == region,]
        
        ExpressionGenes <- rbind(ExpressionGenes, data.frame(Region = rep(region, time = 2),
                                                             `Mod.TF` = rep(mod, times = 2),
                                                             Proportion = c(mean(controlGenes3$Proportion), mean(Rgenes3$Proportion)),
                                                             axisGroup = rep(unique(controlGenes3$axisGroup), times = 2),
                                                             Expression = rep(unique(controlGenes3$Expression), times = 2),
                                                             dataToAnalyse = c("Control gene", "R-gene")))
      }
    }
  }
  
  
  ExpressionGenes$Comparison <- paste(ExpressionGenes$dataToAnalyse, ExpressionGenes$Expression, sep = " \n")
  allGenes$Comparison <- paste(allGenes$dataToAnalyse, allGenes$Expression, sep = " \n")
  
  ExpressionGenes <- ExpressionGenes[order(factor(ExpressionGenes$Comparison, levels = c("Control gene \nNo Expression", "R-gene \nNo Expression", 
                                                                                         "Control gene \nLow Expression", "R-gene \nLow Expression"))),]
  
  allGenes <- allGenes[order(factor(allGenes$Comparison, levels = c("Control gene \nNo Expression", "R-gene \nNo Expression", 
                                                                    "Control gene \nLow Expression", "R-gene \nLow Expression"))),]
  
  
  my_comparisons <- list(c("Control gene \nNo Expression", "R-gene \nNo Expression"), 
                         c("Control gene \nLow Expression", "R-gene \nLow Expression"), 
                         c("R-gene \nNo Expression", "R-gene \nLow Expression"))
  
  for (mod in unique(allGenes$`Mod.TF`)) {
    df <- ExpressionGenes[ExpressionGenes$`Mod.TF`==mod,]
    
    comparison_df <- allGenes[allGenes$`Mod.TF`==mod,]
    
    stat.test <- comparison_df %>% group_by(axisGroup) %>% 
      t_test(Proportion ~ Comparison, comparisons = my_comparisons) %>% 
      mutate(y.position = rep(c(0.9, 0.9, 0.97), times = 10))
    
    plot <- ggbarplot(df, x = "Comparison", y="Proportion", ylab = "Average proportion of gene region",
                      color = "black", fill = "Comparison", 
                      palette = c("azure3", "cadetblue", "bisque2", "darksalmon"), title = mod) + 
      stat_pvalue_manual(
        stat.test, 
        label = "p.adj.signif", size = 4,
        tip.length = 0.01, hide.ns = FALSE) +
      coord_cartesian(ylim= c(0,1), clip = "off") +
      theme_bw() +
      
      geom_text(data = ) +
      
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
    
    
    if (normalised == FALSE) {
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\Non-normalised\\", mod, ".png", sep = ""), plot = plot, width = 10, height = 4)  
    } else if (normalised == TRUE) {
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\Normalised\\", mod, ".png", sep = ""), plot = plot, width = 10, height = 4)  
    }
  }
}