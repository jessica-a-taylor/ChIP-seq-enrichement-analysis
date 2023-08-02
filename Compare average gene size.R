library(ggplot2)
library(data.table)
library(grid)
library(readr)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(rstudioapi)

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
allResultsProportions <- data.frame(read_csv(paste("PlantExp data\\allResultsProportions.csv", sep = "")))
allResultsProportions$Size <- rep("size", times = nrow(allResultsProportions))


# Small genes
small_Rgenes <- geneWidth[which(geneWidth$GeneSet=="R-gene" &
                                  geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25)),]
small_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                        geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25)),]

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
                                  geneWidth$GeneWidth > quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25) &
                                  geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]
medium_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                        geneWidth$GeneWidth > quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .25) &
                                        geneWidth$GeneWidth <= quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]

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
                                   geneWidth$GeneWidth > quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]
big_ControlGenes <- geneWidth[which(geneWidth$GeneSet=="Control gene" &
                                         geneWidth$GeneWidth > quantile(geneWidth[which(geneWidth$GeneSet=="R-gene"),"GeneWidth"], probs = .75)),]

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
jobRunScript("Plot enrichment by size.R",  importEnv = TRUE)