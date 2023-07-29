source("Functions\\Overlaps functions.R")

region = c("UpstreamIntergenic", "Promotor1000", "Promotor500",
           "Gene20", "Gene40", "Gene60", "Gene80", "Gene100", 
           "Downstream", "DownstreamIntergenic")

# Determine the proportion of each gene region overlapping with a significant peak.
frequenciesFunction <- function (allResultsProportions, geneFrequency, geneCount) {

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
                              grepl(level, geneFrequency$Comparison)==TRUE &
                              grepl("Control", geneFrequency$Comparison)==TRUE), "Count"] <- signif((nrow(control_df[which(control_df$Proportion > 0),])/sum(geneCount[which(grepl(level, geneCount$GeneSet)==TRUE &
                                                                                                                                                                  grepl("control", geneCount$GeneSet)==TRUE),"GeneCount"]))*100, digits = 3)
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              grepl(level, geneFrequency$Comparison)==TRUE &
                              grepl("Control", geneFrequency$Comparison)==TRUE), "Enrichment.mean"] <- mean(control_df$Proportion)
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              grepl(level, geneFrequency$Comparison)==TRUE &
                              grepl("Control", geneFrequency$Comparison)==TRUE), "Enrichment.variance"] <- sd(control_df$Proportion)
        
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              grepl(level, geneFrequency$Comparison)==TRUE &
                              grepl("R-gene", geneFrequency$Comparison)==TRUE), "Count"] <- signif((nrow(NLR_df[which(NLR_df$Proportion > 0),])/sum(geneCount[which(grepl(level, geneCount$GeneSet)==TRUE &
                                                                                                                                                                                    grepl("NLR", geneCount$GeneSet)==TRUE),"GeneCount"]))*100, digits = 3)
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              grepl(level, geneFrequency$Comparison)==TRUE &
                              grepl("R-gene", geneFrequency$Comparison)==TRUE), "Enrichment.mean"] <- mean(NLR_df$Proportion)
        
        geneFrequency[which(geneFrequency$Region==r & 
                              geneFrequency$Mod.TF==mod & 
                              grepl(level, geneFrequency$Comparison)==TRUE &
                              grepl("R-gene", geneFrequency$Comparison)==TRUE), "Enrichment.variance"] <- sd(NLR_df$Proportion)
        
      }
    }
  }
  return(geneFrequency)
}
