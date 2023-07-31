# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(genomicData, NLR_genes) {

  # Create a list of files in the folder for a particular tissue type.
  filenames <- list.files(path = "PlantExp data\\Control\\" ,pattern="*.tsv")
  expressionData <- data.frame()
  
  # Merge the data from all files.
  for (file in filenames) {
    expressionData <- rbind(expressionData, as.data.frame(read_tsv(paste("PlantExp data\\Control\\", file, sep = ""), show_col_types = FALSE)))
  }

  # Sample 1000 random genes.
  control_genes <- genomicData[c(sample(nrow(genomicData), 1000)),]
    
  control_data <- expressionData[c(which(expressionData$geneId %in% control_genes$Gene)),]
  control_data <- control_data[,-c(2:4)]
  control_data <- control_data[order(control_data$geneId),]
  control_data$GeneSet <- rep("Control gene", times = nrow(control_data))
  
  # Extract the R-genes.
  NLR_data <- expressionData[c(which(expressionData$geneId %in% NLR_genes$Gene)),]
  NLR_data <- NLR_data[,-c(2:4)]
  NLR_data <- NLR_data[order(NLR_data$geneId),]
  NLR_data$GeneSet <- rep("R-gene", times = nrow(NLR_data))
  
  # Merge the two datasets.
  allGenes <- control_data
  allGenes <- rbind(allGenes, NLR_data)
  
  PlantExpData <- c()
  
  # For each gene, calculate the mean expression across experiments.
  for (gene in unique(allGenes$geneId)) {
    df <- allGenes[allGenes$geneId==gene,]
    
    
    PlantExpData <- rbind(PlantExpData, data.frame(Gene = gene,
                                                   TPM = mean(df$TPM),
                                                   GeneSet = df[1,"GeneSet"]))
  }
  
  # Add a column to PlantExpData with the expression level of each gene.
  expressionLevel <- c()
    
  for (row in 1:nrow(PlantExpData)) {
    if (0 <= PlantExpData[row, "TPM"] & PlantExpData[row, "TPM"] <= 1) {
      expressionLevel <- append(expressionLevel, "No Expression")
    }
    else if (1 < PlantExpData[row, "TPM"] & PlantExpData[row, "TPM"] <= 5) {
      expressionLevel <- append(expressionLevel, "Low Expression")
    }
    else if (5 < PlantExpData[row, "TPM"] & PlantExpData[row, "TPM"] <= 10) {
      expressionLevel <- append(expressionLevel, "Intermediate Expression")
    }
    else if (PlantExpData[row, "TPM"] > 10) {
      expressionLevel <- append(expressionLevel, "High Expression")
    }
  }
  
  PlantExpData <- cbind(PlantExpData, data.frame(ExpressionLevel = expressionLevel))
  
  # Add genomic data to the `PlantExpData` dataframe.
  infoNeeded <- data.frame()
  
  for (row in 1:nrow(PlantExpData)) {
    if (PlantExpData[row, "GeneSet"] == "Control gene") {
      infoNeeded <- rbind(infoNeeded, data.frame(control_genes[which(control_genes$Gene == PlantExpData[row, "Gene"]),]))
    }
    else if (PlantExpData[row, "GeneSet"] == "R-gene") {
      infoNeeded <- rbind(infoNeeded, data.frame(NLR_genes[which(NLR_genes$Gene == PlantExpData[row, "Gene"]),]))
    }
  }
  
  PlantExpData <- cbind(PlantExpData, infoNeeded)
  
  # Sort control and R-genes into a hash based on expression level.
  sampleGenes <- hash()
  
  for (level in unique(PlantExpData$ExpressionLevel)) {
    for (set in c("R-gene", "Control gene")) {
      
      # Randomly sample the same number of control genes as there are R-genes with a particular expression level. 
      if (set == "R-gene") {
        sampleGenes[[paste(set, level, sep = " ")]] <- PlantExpData[which(PlantExpData$Expression == level & PlantExpData$GeneSet == set),]
      }
      else if (set == "Control gene") {
        sampleGenes[[paste(set, level, sep = " ")]] <- PlantExpData[which(PlantExpData$Expression == level & 
                                                                            PlantExpData$GeneSet == set),][c(sample(nrow(PlantExpData[which(PlantExpData$Expression == level & 
                                                                                                  PlantExpData$GeneSet == set),]), nrow(sampleGenes[[paste("R-gene", level, sep = " ")]]))),]
      }
    }
  }
  
  return(sampleGenes)
}
