library(hash)
library(sets)

testData <- data.frame()

source("Get overlapping peaks.R")

allPeaks <- hash()

allModifications <- unique(nextflowOutput$Mod.TF)
allRegions <- names(geneRegions)

i <- 0

totalIterations <- length(geneSet$Gene)

for (gene in geneSet$Gene) {
  
  # Monitor progression - print the fraction of interactions complete.
  i = i + 1
  print(paste(gene, " (", i, "/", totalIterations, ")"))
  
  # Filter for each gene.
  df <- geneSet[geneSet$Gene == gene,]
  
  # Filter for peaks overlapping the gene, promotor and surrounding intergenic regions.
  peaksPerGene <- nextflowOutput[which(nextflowOutput[, "start"] > df[1, "start"] - 5000 &
                                         nextflowOutput[, "end"] < df[1, "end"] + 5000 &
                                         as.numeric(nextflowOutput[, "seqnames"]) == as.numeric(df[1, "seqnames"])),]
  
  # Function to get the peaks for each modification overlapping each gene region.
  allPeaks <- overlappingPeaksFunction(geneRegions, allModifications, allRegions, peaksPerGene, allPeaks, gene)
  
  
  for (mod in allModifications) {
    for (region in allRegions) {
      testData <- rbind(testData, allPeaks[[gene]][[region]][[mod]])
    }
  }
}


write.csv(testData, file = "NewMergedPeaks.csv")