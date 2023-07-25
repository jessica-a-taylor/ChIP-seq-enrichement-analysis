# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(dataToUse, geneRegions) {
  grouping <- c(seq(from = -60, to = -20, by = 20), seq(from = 20, to = 140, by = 20))
  axisGroup <- c()
  
  for (c in 1:length(geneRegions)) {
    axisGroup <- append(axisGroup, rep(grouping[c], times = nrow(dataToUse[dataToUse$Region == geneRegions[c],])))
  }
  dataToUse <- cbind(dataToUse, axisGroup)
  
  return(dataToUse) 
}