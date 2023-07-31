source("Functions\\Overlaps functions.R")

region = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
           "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
           "Downstream", "DownstreamIntergenic")

# Determine the proportion of each gene region overlapping with a significant peak.
frequenciesFunction <- function (allResultsProportions, geneFrequency, geneCount) {

  for (r in unique(allResultsProportions$Region)) {
    for (mod in unique(allResultsProportions$Mod.TF)) {
      for (set in unique(allResultsProportions$GeneSet)) {
        df <- allResultsProportions[which(allResultsProportions$Region==r &
                                            allResultsProportions$Mod.TF==mod &
                                            allResultsProportions$GeneSet==set),]
        

        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$GeneSet==set), "Count"] <- signif((nrow(df[which(df$Proportion > 0),])/sum(geneCount[which(geneCount$GeneSet==set), "GeneCount"]))*100, digits = 2)
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$GeneSet==set), "Enrichment.mean"] <- mean(df$Proportion)
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              geneFrequency$GeneSet==set), "Enrichment.variance"] <- sd(df$Proportion)
      }
    }
  }
  return(geneFrequency)
}
