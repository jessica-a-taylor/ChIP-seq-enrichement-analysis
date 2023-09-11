# Function of determining the coordinates of the intergenic regions between the current and adjacent gene, 
# accounting for their orientation on the DNA strands.

intergenicCoordinatesFunction <- function(dataToUse, genomicData, region) {
  source("Functions\\Get range - merge gene coordinates.R")
  
  genomicData[which(genomicData$strand=="+"), "strand"] <- "positive"
  genomicData[which(genomicData$strand=="-"), "strand"] <- "negative"
  
  # Variations of the functions depending on gene orientations.
  getIntergenicCoordinates <- list(UpstreamIntergenic = 
                                     list(current_positive =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene - 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[adjacentGene, "end"] + 201
                                                     df[which(df$Gene==gene), "end"] <- genomicData[currentGene, "start"] - 1001
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)}),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[adjacentGene, "end"] + 1001
                                                     df[which(df$Gene==gene), "end"] <- genomicData[currentGene, "start"] - 1001
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)})),
                                          
                                          current_negative =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene + 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[currentGene, "end"] + 1001
                                                     df[which(df$Gene==gene), "end"] <- genomicData[adjacentGene, "start"] - 1001
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)}),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[currentGene, "end"] + 1001
                                                     df[which(df$Gene==gene), "end"] <- genomicData[adjacentGene, "start"] - 201
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)}))),
                                   DownstreamIntergenic = 
                                     list(current_positive =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene + 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[currentGene, "end"] + 201
                                                     df[which(df$Gene==gene), "end"] <- genomicData[adjacentGene, "start"] - 1001
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)}),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[currentGene, "end"] + 201
                                                     df[which(df$Gene==gene), "end"] <- genomicData[adjacentGene, "start"] - 201
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)})),
                                          
                                          current_negative =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene - 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[adjacentGene, "end"] + 201
                                                     df[which(df$Gene==gene), "end"] <- genomicData[currentGene, "start"] - 201
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)}),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, df, genomicData) {
                                                     df[which(df$Gene==gene), "start"] <- genomicData[adjacentGene, "end"] + 1001
                                                     df[which(df$Gene==gene), "end"] <- genomicData[currentGene, "start"] - 201
                                                     df[which(df$Gene==gene), "width"] <- df[which(df$Gene==gene), "end"] - df[which(df$Gene==gene), "start"]
                                                     return(df)})))) 
  
  # Determine coordinates of intergenic regions.
  if (nrow(dataToUse) >= 1) {
    df <- genomicData[which(genomicData$Gene %in% dataToUse$Gene),]
    
    df[which(df$strand=="+"), "strand"] <- "positive"
    df[which(df$strand=="-"), "strand"] <- "negative"
  
    for (gene in df$Gene) {
      
      currentGene <- which(genomicData$Gene==gene)
      currentStrand <- paste("current_", genomicData[currentGene,"strand"], sep = "")
      
      adjacentGene <- getIntergenicCoordinates[[region]][[currentStrand]]$getAdjacentGene(currentGene)
      adjacentStrand <- paste("adjacent_", genomicData[adjacentGene,"strand"], sep = "")
      
      df <- getIntergenicCoordinates[[region]][[currentStrand]][[adjacentStrand]]$getCoordinates(currentGene, adjacentGene, df, genomicData)
    }
    df <- df[-which(df$width <= 0),]
      
    df$ranges <- mergeCoordinates(df)
  }
  return(df)
}

# Function for determining the coordinates of the promotor regions 500 bp and 1000 bp upstream of the TSS.
promotorCoordinatesFunction <- function(dataToUse, region) {
  source("Functions\\Get range - merge gene coordinates.R")
  
  df <- dataToUse
  
  df[which(df$strand=="+"), "strand"] <- "positive"
  df[which(df$strand=="-"), "strand"] <- "negative"
  
  # Variations of the functions depending on gene orientations.
  getPromotorCoordinates <- list(Promotor500 =
                                   list(current_positive = 
                                          list(getCoordinates = function(df, dataToUse, row) {
                                            df[row,"start"] <- dataToUse[row,"start"]-500
                                            df[row,"end"] <- dataToUse[row,"start"]
                                            df[row, "width"] <- df[row, "end"] - df[row, "start"]
                                            return(df)}),
                                     
                                     current_negative = 
                                          list(getCoordinates = function(df, dataToUse, row) {
                                            df[row,"start"] <- dataToUse[row,"end"]
                                            df[row,"end"] <- dataToUse[row,"end"]+500
                                            df[row, "width"] <- df[row, "end"] - df[row, "start"]
                                            return(df)})),
                                 
                                 Promotor1000 =
                                   list(current_positive = 
                                          list(getCoordinates = function(df, dataToUse, row) {
                                            df[row,"start"] <- dataToUse[row,"start"]-1000
                                            df[row,"end"] <- dataToUse[row,"start"]
                                            df[row, "width"] <- df[row, "end"] - df[row, "start"]
                                            return(df)}),
                                     
                                        current_negative = 
                                          list(getCoordinates = function(df, dataToUse, row) {
                                            df[row,"start"] <- dataToUse[row,"end"]
                                            df[row,"end"] <- dataToUse[row,"end"]+1000
                                            df[row, "width"] <- df[row, "end"] - df[row, "start"]
                                            return(df)})))
  
  
  if (nrow(df) >= 1) {
    for (row in 1:nrow(df)) {
      currentStrand <- paste("current_", df[row,"strand"], sep = "")
      
      df <- getPromotorCoordinates[[region]][[currentStrand]]$getCoordinates(df, dataToUse, row)
    }
    df$ranges <- mergeCoordinates(df)
  }
  return(df)
}

# Function for determining the coordinates of the 200 bp downstream regions.
downstreamCoordinatesFunction <- function(dataToUse) {
  source("Functions\\Get range - merge gene coordinates.R")
  
  df <- dataToUse
  df[which(df$strand=="+"), "strand"] <- "positive"
  df[which(df$strand=="-"), "strand"] <- "negative"
  
  # Variations of the functions depending on gene orientations.
  getDownstreamCoordinates <- list(current_positive = 
                                   list(getCoordinates = function(df, dataToUse, row) {
                                     df[row,"start"] <- dataToUse[row,"end"]
                                     df[row,"end"] <- dataToUse[row,"end"] + 200
                                     df[row, "width"] <- df[row, "end"] - df[row, "start"]
                                     return(df)}),
                                 
                                 current_negative = 
                                   list(getCoordinates = function(df, dataToUse, row) {
                                     df[row,"start"] <- dataToUse[row,"start"] - 200
                                     df[row,"end"] <- dataToUse[row,"start"]
                                     df[row, "width"] <- df[row, "end"] - df[row, "start"]
                                     return(df)}))
  
  
  if (nrow(df) >= 1) {
    for (row in 1:nrow(df)) {
      currentStrand <- paste("current_", df[row,"strand"], sep = "")
      
      df <- getDownstreamCoordinates[[currentStrand]]$getCoordinates(df, dataToUse, row)
    }
    df$ranges <- mergeCoordinates(df)
  }
  return(df)
}


# Function for determining the coordinates of the gene body in 20% intervals.
genebodyCoordinatesFunction <- function(dataToUse, geneRegions) {
  source("Functions\\Get range - merge gene coordinates.R")
  
  dataToUse[which(dataToUse$strand=="+"), "strand"] <- "positive"
  dataToUse[which(dataToUse$strand=="-"), "strand"] <- "negative"
  
  for (region in c("Gene20", "Gene40", "Gene60", "Gene80", "Gene100")) {
    geneRegions[[region]] <- dataToUse
  }
  
  if (nrow(dataToUse) >= 1) {
    for (row in 1:nrow(dataToUse)) {
      if (dataToUse[row, "strand"]=="positive") {
        geneRegions[["Gene20"]][row, c(4:6)] <- c(dataToUse[row,"start"], 
                                                  dataToUse[row,"start"] + dataToUse[row,"width"]*0.2, 
                                                  (dataToUse[row,"start"] + dataToUse[row,"width"]*0.2) - dataToUse[row,"start"])
        
        geneRegions[["Gene40"]][row, c(4:6)] <- c(geneRegions[["Gene20"]][row,"end"] + 1, 
                                                  geneRegions[["Gene20"]][row,"end"] + 1 + dataToUse[row,"width"]*0.2, 
                                                  (geneRegions[["Gene20"]][row,"end"] + 1 + dataToUse[row,"width"]*0.2) - (geneRegions[["Gene20"]][row,"end"] + 1))
        
        geneRegions[["Gene60"]][row, c(4:6)] <- c(geneRegions[["Gene40"]][row,"end"] + 1, 
                                                  geneRegions[["Gene40"]][row,"end"] + 1 + dataToUse[row,"width"]*0.2, 
                                                  (geneRegions[["Gene40"]][row,"end"] + 1 + dataToUse[row,"width"]*0.2) - (geneRegions[["Gene40"]][row,"end"] + 1))
        
        geneRegions[["Gene80"]][row, c(4:6)] <- c(geneRegions[["Gene60"]][row,"end"] + 1, 
                                                  geneRegions[["Gene60"]][row,"end"] + 1 + dataToUse[row,"width"]*0.2, 
                                                  (geneRegions[["Gene60"]][row,"end"] + 1 + dataToUse[row,"width"]*0.2) - (geneRegions[["Gene60"]][row,"end"] + 1))
        
        geneRegions[["Gene100"]][row, c(4:6)] <- c(geneRegions[["Gene80"]][row,"end"] + 1, 
                                                   dataToUse[row,"end"], 
                                                   dataToUse[row,"end"] - (geneRegions[["Gene80"]][row,"end"] + 1))
      }
      else if (dataToUse[row, "strand"]=="negative"){
        geneRegions[["Gene20"]][row, c(4:6)] <- c(dataToUse[row,"end"] - dataToUse[row,"width"]*0.2, 
                                                  dataToUse[row,"end"], 
                                                  dataToUse[row,"end"] - (dataToUse[row,"end"] - dataToUse[row,"width"]*0.2))
        
        geneRegions[["Gene40"]][row, c(4:6)] <- c(geneRegions[["Gene20"]][row,"start"] - dataToUse[row,"width"]*0.2, 
                                                  geneRegions[["Gene20"]][row,"start"] - 1, 
                                                  (geneRegions[["Gene20"]][row,"start"] - 1) - (geneRegions[["Gene20"]][row,"start"] - dataToUse[row,"width"]*0.2))
        
        geneRegions[["Gene60"]][row, c(4:6)] <- c(geneRegions[["Gene40"]][row,"start"] - dataToUse[row,"width"]*0.2, 
                                                  geneRegions[["Gene40"]][row,"start"] - 1, 
                                                  (geneRegions[["Gene40"]][row,"start"] - 1) - (geneRegions[["Gene40"]][row,"start"] - dataToUse[row,"width"]*0.2))
        
        geneRegions[["Gene80"]][row, c(4:6)] <- c(geneRegions[["Gene60"]][row,"start"] - dataToUse[row,"width"]*0.2, 
                                                  geneRegions[["Gene60"]][row,"start"] - 1, 
                                                  (geneRegions[["Gene60"]][row,"start"] - 1) - (geneRegions[["Gene60"]][row,"start"] - dataToUse[row,"width"]*0.2))
        
        geneRegions[["Gene100"]][row, c(4:6)] <- c(dataToUse[row,"start"],
                                                   geneRegions[["Gene80"]][row,"start"] - 1, 
                                                   (geneRegions[["Gene80"]][row,"start"] - 1) - dataToUse[row,"start"])
      }
    }
  }
  return(geneRegions)
}