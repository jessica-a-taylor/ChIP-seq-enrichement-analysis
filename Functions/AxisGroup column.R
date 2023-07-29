# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(dataToUse, geneRegions) {
  grouping <- c(seq(from = -60, to = -20, by = 20), seq(from = 20, to = 100, by = 20), seq(from = 140, to = 160, by = 20))
  axisGroup <- c("UpstreamIntergenic", "Promotor1000", "Promotor500", "Gene20", "Gene40", "Gene60",
                 "Gene80", "Gene100", "Downstream", "DownstreamIntergenic")
  
  # Add empty 'axisGroup' column to 'dataToUse.
  dataToUse <- cbind(dataToUse, data.frame(axisGroup = rep(NA, times = nrow(dataToUse))))
  
  for (c in 1:length(names(geneRegions))) {
    dataToUse[which(dataToUse$Region == axisGroup[c]),"axisGroup"] <- grouping[c]
  }
  return(dataToUse) 
}
