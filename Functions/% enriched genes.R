source("Functions\\Overlaps functions.R")

region = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
           "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
           "Downstream", "DownstreamIntergenic")

# Determine the proportion of each gene region overlapping with a significant peak.
frequenciesFunction <- function (allResultsProportions, geneFrequency) {

  for (r in unique(allResultsProportions$Region)) {
    for (mod in unique(allResultsProportions$Mod.TF)) {
      for (level in unique(allResultsProportions$Expression)) {
        df <- allResultsProportions[which(allResultsProportions$Region==r &
                                            allResultsProportions$Mod.TF==mod &
                                            allResultsProportions$Expression==level),]
        
        control_df <- df[which(grepl("control", df$dataToAnalyse) == TRUE),]
        NLR_df <- df[which(grepl("NLR", df$dataToAnalyse) == TRUE),]

        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$Expression==level &
                              geneFrequency$GeneSet=="Control gene"), "Count"] <- nrow(control_df[which(control_df$Proportion > 0),])
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$Expression==level &
                              geneFrequency$GeneSet=="Control gene"), "Enrichment.variance"] <- sd(control_df$Proportion)
        
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$Expression==level &
                              geneFrequency$GeneSet=="R-gene"), "Count"] <- nrow(NLR_df[which(NLR_df$Proportion > 0),])
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$Expression==level &
                              geneFrequency$GeneSet=="R-gene"), "Enrichment.variance"] <- sd(NLR_df$Proportion)
        
      }
    }
  }
  return(geneFrequency)
}
