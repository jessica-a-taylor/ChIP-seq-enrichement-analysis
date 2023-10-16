# Function to add a column to the dataframe with the numbers on the x axis that will correspond with each gene region.
geneRegionAxisLocations <- function(proportion, geneRegions, wb) {
  grouping <- c(seq(from = -60, to = -20, by = 20), seq(from = 20, to = 100, by = 20), seq(from = 140, to = 160, by = 20))
  axisGroup <- c("UpstreamIntergenic", "Promotor1000", "Promotor500", "Gene20", "Gene40", "Gene60",
                 "Gene80", "Gene100", "Downstream", "DownstreamIntergenic")
  
  allModifications <- unique(nextflowOutput$Mod.TF)
  allRegions <- names(geneRegions)
  
  for (mod in allModifications) {
    addWorksheet(wb, sheetName = mod)
    wb_data <- data.frame()
   
    for (gene in names(proportion[[mod]])) {
      for (region in allRegions) {
        if (nrow(proportion[[mod]][[gene]][[region]]) >= 1) {
          
          proportion[[mod]][[gene]][[region]]$axisGroup <- grouping[which(axisGroup == region)]
          wb_data <- rbind(wb_data, proportion[[mod]][[gene]][[region]])
        }
      }
    } 
    writeData(wb, sheet = mod, wb_data)
  }
  return(wb) 
}
