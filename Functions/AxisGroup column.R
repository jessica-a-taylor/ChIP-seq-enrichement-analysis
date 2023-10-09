# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(proportion, geneRegions) {
  grouping <- c(seq(from = -60, to = -20, by = 20), seq(from = 20, to = 100, by = 20), seq(from = 140, to = 160, by = 20))
  axisGroup <- c("UpstreamIntergenic", "Promotor1000", "Promotor500", "Gene20", "Gene40", "Gene60",
                 "Gene80", "Gene100", "Downstream", "DownstreamIntergenic")
  
  allModifications <- unique(nextflowOutput$Mod.TF)
  allRegions <- names(geneRegions)
  
  for (gene in names(proportion)) {
    for (region in allRegions) {
      for (mod in allModifications) {
        if (nrow(proportion[[gene]][[region]][[mod]]) >= 1) {
          
         for (c in 1:length(allRegions)) {
            proportion[[gene]][[region]][[mod]]$axisGroup <- grouping[which(axisGroup == region)]
          }
        }
      }
    }
  }
  return(proportion) 
}
