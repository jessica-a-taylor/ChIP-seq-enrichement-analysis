# Get the coordinates of each gene region.
source("Functions\\getCoordinates.R")

geneRegions2 <- hash()

# Determine the coordinates of the intergenic regions.
for (region in c("UpstreamIntergenic", "DownstreamIntergenic")) {
  geneRegions2[[region]] <- intergenicCoordinatesFunction(dataToUse, genomicData, region)
}

# Determine the coordinates of the promotor regions.
for (region in c("Promotor500", "Promotor1000")) {
  geneRegions2[[region]] <- promotorCoordinatesFunction(dataToUse, region)
}

# Determine the coordinates of the 200 bp downstream region.
geneRegions2[["Downstream"]] <- downstreamCoordinatesFunction(dataToUse)

# Determine the coordinates of the gene body in 20% intervals.
geneRegions2 <- genebodyCoordinatesFunction(dataToUse, geneRegions2)