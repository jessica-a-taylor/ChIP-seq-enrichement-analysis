library(ggplot2)
library(data.table)
library(grid)
library(readr)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(rstudioapi)

normalised <- FALSE

# Compare average gene size between R-genes and control genes.
geneWidth <- data.frame()

for (set in names(sampleGenesPlantExp)) {
  geneWidth <- rbind(geneWidth, data.frame(Gene = sampleGenesPlantExp[[set]]$Gene,
                                           GeneSet = sampleGenesPlantExp[[set]]$GeneSet,
                                           GeneWidth = sampleGenesPlantExp[[set]]$width/1000))
  
}


plot <- ggplot(geneWidth, aes(x = GeneSet, y = GeneWidth)) +
  geom_boxplot() + labs(x = "Gene set", y = "Gene width (kb)") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "R-gene") + theme_bw() 

ggsave("Graphs\\Gene width comparison.png", plot = plot, width = 8, height = 4)  


# Determine whether there is still a significant difference between R-genes and control genes of each size category.
allResultsProportions <- data.frame(read_csv(paste("PlantExp data\\Non-normalised\\allResultsProportions.csv", sep = "")))
allResultsProportions$Size <- rep("size", times = nrow(allResultsProportions))

# Very small genes
verySmall_Rgenes <- geneWidth[which(geneWidth$GeneSet=="R-gene" & 
                                  geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .10)),]
verySmall_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                        geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .10)),]

if (nrow(verySmall_ControlGenes) >= nrow(verySmall_Rgenes)) {
  verySmall_ControlGenes <- verySmall_ControlGenes[c(sample(nrow(verySmall_ControlGenes), nrow(verySmall_Rgenes))),]
} else   verySmall_ControlGenes <- verySmall_ControlGenes

verySmallGenes <- verySmall_Rgenes
verySmallGenes <- rbind(verySmallGenes, verySmall_ControlGenes)

plot <- ggplot(verySmallGenes, aes(x = GeneSet, y = GeneWidth)) +
  geom_boxplot() + labs(x = paste("Gene set (n = ", nrow(verySmall_Rgenes), ")", sep = ""), y = "Gene width (kb)") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "R-gene") + theme_bw()

ggsave("Graphs\\Very small gene width comparison.png", plot = plot, width = 8, height = 4)  

# Add gene size to 'allResultsProportions'.
for (row in 1:nrow(verySmallGenes)) {
  allResultsProportions[which(allResultsProportions$Gene == verySmallGenes[row,"Gene"]),"Size"] <- rep("Very small", times = nrow(allResultsProportions[which(allResultsProportions$Gene == verySmallGenes[row,"Gene"]),]))
}


# Small genes
small_Rgenes <- geneWidth[which(geneWidth$GeneSet=="R-gene" &
                                  geneWidth$GeneWidth >= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .10) &
                                  geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25)),]
small_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                        geneWidth$GeneWidth >= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .10) &
                                        geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25)),]

if (nrow(small_ControlGenes) >= nrow(small_Rgenes)) {
  small_ControlGenes <- small_ControlGenes[c(sample(nrow(small_ControlGenes), nrow(small_Rgenes))),]
} else   small_ControlGenes <- small_ControlGenes

smallGenes <- small_Rgenes
smallGenes <- rbind(smallGenes, small_ControlGenes)

plot <- ggplot(smallGenes, aes(x = GeneSet, y = GeneWidth)) +
  geom_boxplot() + labs(x = paste("Gene set (n = ", nrow(small_Rgenes), ")", sep = ""), y = "Gene width (kb)") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "R-gene") + theme_bw()

ggsave("Graphs\\Small gene width comparison.png", plot = plot, width = 8, height = 4) 

# Add gene size to 'allResultsProportions'.
for (row in 1:nrow(smallGenes)) {
  allResultsProportions[which(allResultsProportions$Gene == smallGenes[row,"Gene"]),"Size"] <- rep("Small", times = nrow(allResultsProportions[which(allResultsProportions$Gene == smallGenes[row,"Gene"]),]))
}


# Medium genes
medium_Rgenes <- geneWidth[which(geneWidth$GeneSet=="R-gene" &
                                  geneWidth$GeneWidth >= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25) &
                                  geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]
medium_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                        geneWidth$GeneWidth >= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25) &
                                        geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]

if (nrow(medium_ControlGenes) >= nrow(medium_Rgenes)) {
  medium_ControlGenes <- medium_ControlGenes[c(sample(nrow(medium_ControlGenes), nrow(medium_Rgenes))),]
} else   medium_ControlGenes <- medium_ControlGenes

mediumGenes <- medium_Rgenes
mediumGenes <- rbind(mediumGenes, medium_ControlGenes)

plot <- ggplot(mediumGenes, aes(x = GeneSet, y = GeneWidth)) +
  geom_boxplot() + labs(x = paste("Gene set (n = ", nrow(medium_Rgenes), ")", sep = ""), y = "Gene width (kb)") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "R-gene") + theme_bw()

ggsave("Graphs\\Medium gene width comparison.png", plot = plot, width = 8, height = 4) 

# Add gene size to 'allResultsProportions'.
for (row in 1:nrow(mediumGenes)) {
  allResultsProportions[which(allResultsProportions$Gene == mediumGenes[row,"Gene"]),"Size"] <- rep("Medium", times = nrow(allResultsProportions[which(allResultsProportions$Gene == mediumGenes[row,"Gene"]),]))
}


# Big genes
big_Rgenes <- geneWidth[which(geneWidth$GeneSet=="R-gene" &
                                   geneWidth$GeneWidth >= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]
big_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                         geneWidth$GeneWidth >= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]

if (nrow(big_ControlGenes) >= nrow(big_Rgenes)) {
  big_ControlGenes <- big_ControlGenes[c(sample(nrow(big_ControlGenes), nrow(big_Rgenes))),]
} else   big_ControlGenes <- big_ControlGenes

bigGenes <- big_Rgenes
bigGenes <- rbind(bigGenes, big_ControlGenes)

plot <- ggplot(bigGenes, aes(x = GeneSet, y = GeneWidth)) +
  geom_boxplot() + labs(x = paste("Gene set (n = ", nrow(big_Rgenes), ")", sep = ""), y = "Gene width (kb)") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "R-gene") + theme_bw()

ggsave("Graphs\\Big gene width comparison.png", plot = plot, width = 8, height = 4) 

# Add gene size to 'allResultsProportions'.
for (row in 1:nrow(bigGenes)) {
  allResultsProportions[which(allResultsProportions$Gene == bigGenes[row,"Gene"]),"Size"] <- rep("Large", times = nrow(allResultsProportions[which(allResultsProportions$Gene == bigGenes[row,"Gene"]),]))
}

# Plot bargraph.
jobRunScript("Functions\\Enrichment plot function.R",  importEnv = TRUE)


# Repeat the enrichment analysis for genes in each size category.
allResultsProportions <- data.frame(read_csv(paste(analysis, "\\Non-normalised\\allResultsProportions.csv", sep = "")))

axisText <- c("Intergenic", "Promotor \n(1kb)", "Promotor \n(500bp)", "TSS", "20%",
              "40%", "60%", "80%", "100%", "TTS", "Downstream \n(200bp)", "Intergenic")



# Enrichment analysis for small genes with widths <= the average width of the R-genes (3.868219 kb).
small_Rgenes <- geneWidth[which(grepl("R-genes", geneWidth$GeneSet)==TRUE & geneWidth$GeneWidth <= 3.868219),"Gene"]
small_ControlGenes <- geneWidth[which(grepl("control", geneWidth$GeneSet)==TRUE & geneWidth$GeneWidth <= 3.868219),"Gene"]


allGenes <- allResultsProportions[which(allResultsProportions$Expression == "No Expression" |
                                          allResultsProportions$Expression == "Low Expression"),]
 
controlGenes <- allGenes[grepl("control", allGenes$dataToAnalyse),]
controlGenes <- controlGenes[which(controlGenes$Gene %in% small_ControlGenes),]
controlGenes$dataToAnalyse <- rep("Control gene", times = nrow(controlGenes))
  
Rgenes <- allGenes[grepl("NLR", allGenes$dataToAnalyse),]
Rgenes <- Rgenes[which(Rgenes$Gene %in% small_Rgenes),]
Rgenes$dataToAnalyse <- rep("R-gene", times = nrow(Rgenes))
  
allGenes <- controlGenes
allGenes <- rbind(allGenes, Rgenes)
  
ExpressionGenes <- data.frame()
  
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
  
for (mod in unique(allGenes$Mod.TF)) {
  df <- ExpressionGenes[ExpressionGenes$`Mod.TF`==mod,]
  
  comparison_df <- allGenes[allGenes$`Mod.TF`==mod,]
  
  stat.test <- comparison_df %>% group_by(axisGroup) %>% 
    t_test(Proportion ~ Comparison, comparisons = my_comparisons) %>% 
    mutate(y.position = rep(c(0.92, 0.92, 0.98), times = 10))
  
  plot <- ggbarplot(df, x = "Comparison", y="Proportion", ylab = "Average enrichment",
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
  
  plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                font.ytickslab = 8)
  
  ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", mod, "-small-genes.png", sep = ""), plot = plot, width = 10, height = 4)  
} 


# Enrichment analysis for big genes with widths >= the average width of the R-genes (3.868219 kb).
big_Rgenes <- geneWidth[which(grepl("R-genes", geneWidth$GeneSet)==TRUE & geneWidth$GeneWidth > 3.868219),"Gene"]
big_ControlGenes <- geneWidth[which(grepl("control", geneWidth$GeneSet)==TRUE & geneWidth$GeneWidth > 3.868219),"Gene"]
  
allGenes <- allResultsProportions[which(allResultsProportions$Expression == "No Expression" |
                                          allResultsProportions$Expression == "Low Expression"),]

controlGenes <- allGenes[grepl("control", allGenes$dataToAnalyse),]
controlGenes <- controlGenes[which(controlGenes$Gene %in% big_ControlGenes),]
controlGenes$dataToAnalyse <- rep("Control gene", times = nrow(controlGenes))

Rgenes <- allGenes[grepl("NLR", allGenes$dataToAnalyse),]
Rgenes <- Rgenes[which(Rgenes$Gene %in% big_Rgenes),]
Rgenes$dataToAnalyse <- rep("R-gene", times = nrow(Rgenes))

allGenes <- controlGenes
allGenes <- rbind(allGenes, Rgenes)

ExpressionGenes <- data.frame()

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
                                                             axisGroup = rep(unique(Rgenes3$axisGroup), times = 2),
                                                             Expression = rep(unique(Rgenes3$Expression), times = 2),
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

for (mod in unique(allGenes$Mod.TF)) {
  df <- ExpressionGenes[ExpressionGenes$`Mod.TF`==mod,]
  
  comparison_df <- allGenes[allGenes$`Mod.TF`==mod,]
  
  stat.test <- comparison_df %>% group_by(axisGroup) %>% 
    t_test(Proportion ~ Comparison, comparisons = my_comparisons) %>% 
    mutate(y.position = rep(c(0.92, 0.92, 0.98), times = 10))
  
  plot <- ggbarplot(df, x = "Comparison", y="Proportion", ylab = "Average enrichment",
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
  
  plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                font.ytickslab = 8)
  
  ggsave(paste("Graphs\\Enrichment\\", analysis, "\\", mod, "-big-genes.png", sep = ""), plot = plot, width = 10, height = 4)  
} 