# Function to add a column for expression level.
expressionColumn <- function(proportion, set) {
  
  allModifications <- unique(nextflowOutput$Mod.TF)
  allRegions <- names(geneRegions)
  
  for (gene in names(proportion)) {
    for (region in allRegions) {
      for (mod in allModifications) {
        if (nrow(proportion[[gene]][[region]][[mod]]) >= 1) {
          
          proportion[[gene]][[region]][[mod]] <- cbind(proportion[[gene]][[region]][[mod]], 
                                                       data.frame(GeneSet = rep(set, times = nrow(proportion[[gene]][[region]][[mod]]))))
        }
      }
    }
  }
  return(proportion) 
}