# Get the coordinates of each gene region.
getGeneCoordinates <- function(dataToUse, genomicData) {
  source("Functions\\getCoordinates.R")
  
  geneRegions <- hash()
  
  # Determine the coordinates of the intergenic regions.
  for (region in c("UpstreamIntergenic", "DownstreamIntergenic")) {
    geneRegions[[region]] <- intergenicCoordinatesFunction(dataToUse, genomicData, region)
  }
  
  # Determine the coordinates of the promotor regions.
  for (region in c("Promotor500", "Promotor1000")) {
    geneRegions[[region]] <- promotorCoordinatesFunction(dataToUse, region)
  }
  
  # Determine the coordinates of the 200 bp downstream region.
  geneRegions[["Downstream"]] <- downstreamCoordinatesFunction(dataToUse)
  
  # Determine the coordinates of the gene body in 20% intervals.
  geneRegions <- genebodyCoordinatesFunction(dataToUse, geneRegions)
}
