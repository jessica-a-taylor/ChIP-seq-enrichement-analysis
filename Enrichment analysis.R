# Load required libraries.


# Enrichment analysis based on the occurrence of significant peaks.
setwd("C:/Users/jexy2/OneDrive/ChIP-seq-enrichment-analysis")

# Import ChIP-seq data.
source("Functions\\getChIP_seq_data.R")
nextflowOutput <- getChIP_seq_data()

# Sample 1000 random control genes, then sort R-genes and control genes based on expression level.
# Ensure that the number of R-genes and control genes is the same for a particular expression level.
source("Functions\\PlantExp.R")
sampleGenes <- PlantExp()

# Perform enrichment analysis.
source("Functions\\analysisFunction.R")
analysisFunction(sampleGenes, nextflowOutput)

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

# Generate enrichment plots.
jobRunScript("Plot enrichment.R",  importEnv = TRUE)
jobRunScript("Compare average gene size.R",  importEnv = TRUE)