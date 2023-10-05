# Function of determining the coordinates of the intergenic regions between the current and adjacent gene, 
# accounting for their orientation on the DNA strands.

intergenicCoordinatesFunction <- function(dataToUse, genomicData, region) {
  
  # Variations of the functions depending on gene orientations.
  getIntergenicCoordinates <- list(UpstreamIntergenic = 
                                     list(current_positive =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene - 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[currentGene, "start"] - 1001)-(genomicData[adjacentGene, "end"] + 201) > 0) {
                                                       data[[gene]]$start <- genomicData[adjacentGene, "end"] + 201
                                                       data[[gene]]$end <- genomicData[currentGene, "start"] - 1001
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[currentGene, "start"] - 1001) - (genomicData[adjacentGene, "end"] + 1001) > 0) {
                                                       data[[gene]]$start <- genomicData[adjacentGene, "end"] + 1001
                                                       data[[gene]]$end <- genomicData[currentGene, "start"] - 1001
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     })),
                                          
                                          current_negative =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene + 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[adjacentGene, "start"] - 1001) - (genomicData[currentGene, "end"] + 1001) > 0) {
                                                       data[[gene]]$start <- genomicData[currentGene, "end"] + 1001
                                                       data[[gene]]$end <- genomicData[adjacentGene, "start"] - 1001
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[adjacentGene, "start"] - 201) - (genomicData[currentGene, "end"] + 1001) > 0) {
                                                       data[[gene]]$start <- genomicData[currentGene, "end"] + 1001
                                                       data[[gene]]$end <- genomicData[adjacentGene, "start"] - 201
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data) 
                                                     return(data)
                                                     }))),
                                   
                                   DownstreamIntergenic = 
                                     list(current_positive =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene + 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[adjacentGene, "start"] - 1001) - (genomicData[currentGene, "end"] + 201) > 0) {
                                                       data[[gene]]$start <- genomicData[currentGene, "end"] + 201
                                                       data[[gene]]$end <- genomicData[adjacentGene, "start"] - 1001
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[adjacentGene, "start"] - 201) - (genomicData[currentGene, "end"] + 201) > 0) {
                                                       data[[gene]]$start <- genomicData[currentGene, "end"] + 201
                                                       data[[gene]]$end <- genomicData[adjacentGene, "start"] - 201
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     })),
                                          
                                          current_negative =
                                            list(getAdjacentGene = function(currentGene) {return(currentGene - 1)},
                                                 
                                                 adjacent_positive =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[currentGene, "start"] - 201) - (genomicData[adjacentGene, "end"] + 201) > 0) {
                                                       data[[gene]]$start <- genomicData[adjacentGene, "end"] + 201
                                                       data[[gene]]$end <- genomicData[currentGene, "start"] - 201
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     }),
                                                 
                                                 adjacent_negative =
                                                   list(getCoordinates = function(currentGene, adjacentGene, data, gene, genomicData) {
                                                     if ((genomicData[currentGene, "start"] - 201) - (genomicData[adjacentGene, "end"] + 1001) > 0) {
                                                       data[[gene]]$start <- genomicData[adjacentGene, "end"] + 1001
                                                       data[[gene]]$end <- genomicData[currentGene, "start"] - 201
                                                       data[[gene]]$width <- data[[gene]]$end - data[[gene]]$start
                                                       data[[gene]]$ranges <- paste(data[[gene]]$start, "-", data[[gene]]$end, sep = "")
                                                     } else del(gene, data)
                                                     return(data)
                                                     })))) 
  
  # Determine coordinates of intergenic regions.
  if (nrow(dataToUse) >= 1) {
    
    # Create hash storing the information for each gene.
    geneCoordinates <- hash()
    
    # Change the '+/-' symbols to words so that they can be used as keys in the hash.
    genomicData[which(genomicData$strand=="+"), "strand"] <- "positive"
    genomicData[which(genomicData$strand=="-"), "strand"] <- "negative"
    
    df <- genomicData[which(genomicData$Gene %in% dataToUse$Gene),]
  
    for (gene in df$Gene) {
      
      geneCoordinates[[gene]] <- df[df$Gene==gene,]
      
      currentGene <- which(genomicData$Gene==gene)
      currentStrand <- paste("current_", genomicData[currentGene,"strand"], sep = "")
      
      adjacentGene <- getIntergenicCoordinates[[region]][[currentStrand]]$getAdjacentGene(currentGene)
      adjacentStrand <- paste("adjacent_", genomicData[adjacentGene,"strand"], sep = "")
      
      geneCoordinates <- getIntergenicCoordinates[[region]][[currentStrand]][[adjacentStrand]]$getCoordinates(currentGene, adjacentGene, geneCoordinates, gene, genomicData)
    }
  }
  return(geneCoordinates)
}

# Function for determining the coordinates of the promotor regions 500 bp and 1000 bp upstream of the TSS.
promotorCoordinatesFunction <- function(dataToUse, region) {
  
  # Variations of the functions depending on gene orientations.
  getPromotorCoordinates <- list(Promotor500 =
                                   list(current_positive = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$start-500
                                            data$end <- data$start
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)}),
                                     
                                     current_negative = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$end
                                            data$end <- data$end+500
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)})),
                                 
                                 Promotor1000 =
                                   list(current_positive = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$start-1000
                                            data$end <- data$start
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)}),
                                     
                                        current_negative = 
                                          list(getCoordinates = function(data) {
                                            data$start <- data$end
                                            data$end <- data$end+1000
                                            data$width <- data$end - data$start
                                            data$ranges <- paste(data$start, "-", data$end, sep = "")
                                            return(data)})))
  
  
  # Determine coordinates of promotor regions.
  if (nrow(dataToUse) >= 1) {
    
    # Create hash storing the information for each gene.
    geneCoordinates <- hash()
    
    # Change the '+/-' symbols to words so that they can be used as keys in the hash.
    dataToUse[which(dataToUse$strand=="+"), "strand"] <- "positive"
    dataToUse[which(dataToUse$strand=="-"), "strand"] <- "negative"
    
    for (row in 1:nrow(dataToUse)) {
      geneCoordinates[[dataToUse[row,"Gene"]]] <- dataToUse[row,]
      
      currentStrand <- paste("current_", dataToUse[row,"strand"], sep = "")
      
      geneCoordinates[[dataToUse[row,"Gene"]]] <- getPromotorCoordinates[[region]][[currentStrand]]$getCoordinates(geneCoordinates[[dataToUse[row,"Gene"]]])
    }
  }
  return(geneCoordinates)
}

# Function for determining the coordinates of the 200 bp downstream regions.
downstreamCoordinatesFunction <- function(dataToUse) {
 
  # Variations of the functions depending on gene orientations.
  getDownstreamCoordinates <- list(current_positive = 
                                   list(getCoordinates = function(data) {
                                     data$start <- data$end
                                     data$end <- data$end + 200
                                     data$width <- data$end - data$start
                                     data$ranges <- paste(data$start, "-", data$end, sep = "")
                                     return(data)}),
                                 
                                 current_negative = 
                                   list(getCoordinates = function(data) {
                                     data$start <- data$start - 200
                                     data$end <- data$start
                                     data$width <- data$end - data$start
                                     data$ranges <- paste(data$start, "-", data$end, sep = "")
                                     return(data)}))
  
  
  # Determine coordinates of downstream regions.
  if (nrow(dataToUse) >= 1) {
    
    # Create hash storing the information for each gene.
    geneCoordinates <- hash()
    
    # Change the '+/-' symbols to words so that they can be used as keys in the hash.
    dataToUse[which(dataToUse$strand=="+"), "strand"] <- "positive"
    dataToUse[which(dataToUse$strand=="-"), "strand"] <- "negative"
    
    for (row in 1:nrow(dataToUse)) {
      geneCoordinates[[dataToUse[row,"Gene"]]] <- dataToUse[row,]
      
      currentStrand <- paste("current_", dataToUse[row,"strand"], sep = "")
      
      geneCoordinates[[dataToUse[row,"Gene"]]] <- getDownstreamCoordinates[[currentStrand]]$getCoordinates(geneCoordinates[[dataToUse[row,"Gene"]]])
    }
  }
  return(geneCoordinates)
}


# Function for determining the coordinates of the gene body in 20% intervals.
genebodyCoordinatesFunction <- function(dataToUse, geneRegions) {
  
  if (nrow(dataToUse) >= 1) {
    
    # Create hash storing the information for each gene.
    for (region in c("Gene20", "Gene40", "Gene60", "Gene80", "Gene100")) {
      geneRegions[[region]] <- hash()
    }

    # Change the '+/-' symbols to words so that they can be used as keys in the hash.
    dataToUse[which(dataToUse$strand=="+"), "strand"] <- "positive"
    dataToUse[which(dataToUse$strand=="-"), "strand"] <- "negative"
    
    for (row in 1:nrow(dataToUse)) {
      for (region in c("Gene20", "Gene40", "Gene60", "Gene80", "Gene100")) {
        geneRegions[[region]][[dataToUse[row,"Gene"]]] <- dataToUse[row,]
      }
      
      if (dataToUse[row, "strand"]=="positive") {
        
        geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(dataToUse[row,"start"]), 
                                                                          as.numeric(dataToUse[row,"start"]) + (as.numeric(dataToUse[row,"width"])*0.2),
                                                                          (as.numeric(dataToUse[row,"start"]) + (as.numeric(dataToUse[row,"width"])*0.2) - as.numeric(dataToUse[row,"start"])),
                                                                          paste(dataToUse[row,"start"], "-", dataToUse[row,"start"] + dataToUse[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$end) + 1, 
                                                                          as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$end) + 1 + (as.numeric(dataToUse[row,"width"])*0.2), 
                                                                          (as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$end) + 1 + (as.numeric(dataToUse[row,"width"])*0.2) - (as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$end) + 1)),
                                                                          paste(as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$end) + 1, "-", as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$end) + 1 + dataToUse[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$end) + 1,
                                                                          as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$end) + 1 + (as.numeric(dataToUse[row,"width"])*0.2),
                                                                          (as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$end) + 1 + (as.numeric(dataToUse[row,"width"])*0.2) - (as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$end) + 1)),
                                                                          paste(as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$end) + 1, "-", as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$end) + 1 + dataToUse[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$end) + 1,
                                                                          as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$end) + 1 + (as.numeric(dataToUse[row,"width"])*0.2),
                                                                          (as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$end) + 1 + (as.numeric(dataToUse[row,"width"])*0.2) - (as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$end) + 1)),
                                                                          paste(as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$end) + 1, "-", as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$end) + 1 + dataToUse[row,"width"]*0.2, sep = ""))
        
        geneRegions[["Gene100"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]]$end) + 1, 
                                                                           as.numeric(dataToUse[row,"end"]), 
                                                                           as.numeric(dataToUse[row,"end"]) - (as.numeric(geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]]$end) + 1),
                                                                           paste(as.numeric(geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]]$end) + 1, "-", as.numeric(dataToUse[row,"end"]), sep = ""))
      }
      else if (dataToUse[row, "strand"]=="negative"){
        geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(dataToUse[row,"end"]) - (as.numeric(dataToUse[row,"width"])*0.2),
                                                                          as.numeric(dataToUse[row,"end"]),
                                                                          as.numeric(dataToUse[row,"end"] - (dataToUse[row,"end"] - dataToUse[row,"width"]*0.2)),
                                                                          paste(dataToUse[row,"end"] - dataToUse[row,"width"]*0.2, "-", dataToUse[row,"end"], sep = ""))
        
        geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2), 
                                                                          as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$start) - 1,
                                                                          (as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$start) - 1) - (as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2)),
                                                                          paste(as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2), "-", as.numeric(geneRegions[["Gene20"]][[dataToUse[row,"Gene"]]]$start) - 1, sep = ""))
        
        geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2),
                                                                          as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$start) - 1, 
                                                                          (as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$start) - 1) - (as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2)),
                                                                          paste(as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2), "-", as.numeric(geneRegions[["Gene40"]][[dataToUse[row,"Gene"]]]$start) - 1, sep = ""))
        
        geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2),
                                                                          as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$start) - 1,
                                                                          (as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$start) - 1) - (as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2)),
                                                                          paste(as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$start) - (as.numeric(dataToUse[row,"width"])*0.2), "-", as.numeric(geneRegions[["Gene60"]][[dataToUse[row,"Gene"]]]$start) - 1, sep = ""))
        
        geneRegions[["Gene100"]][[dataToUse[row,"Gene"]]][, c(4:6,8)] <- c(as.numeric(dataToUse[row,"start"]),
                                                                           as.numeric(geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]]$start) - 1,
                                                                           (as.numeric(geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]]$start) - 1) - as.numeric(dataToUse[row,"start"]),
                                                                           paste(as.numeric(dataToUse[row,"start"]), "-", as.numeric(geneRegions[["Gene80"]][[dataToUse[row,"Gene"]]]$start) - 1, sep = ""))
      }
    }
  }
  return(geneRegions)
}