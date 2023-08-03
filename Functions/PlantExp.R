# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(genomicData, NLR_genes) {

  # Import expression data.
  PlantExpData <- as.data.frame(read_csv("PlantExp data\\geneExpression.csv"))
  PlantExpData <- PlantExpData[-c(which(is.na(PlantExpData$TPM))),]
  

  # Extract the R-genes.
  NLR_data <- PlantExpData[c(which(PlantExpData$Gene %in% NLR_genes$Gene)),]
  NLR_data$GeneSet <- rep("R-gene", times = nrow(NLR_data))
  
  
  # Sample 1000 random genes with a gene length and expression distribution equal to that of the R-genes.
  if (normalised == TRUE) {
    control_data <- data.frame()
    candidateGenes <- c()
    for (n in seq(from = .1, to = 1, by = .1)) {
      candidateGenes <- append(candidateGenes, sample(PlantExpData[which(PlantExpData$width > quantile(NLR_data$width, probs = n-.1) &
                                                                           PlantExpData$width <= quantile(NLR_data$width, probs = n-.05) &
                                                                          !(PlantExpData$Gene %in% candidateGenes) &
                                                                           PlantExpData$TPM < max(NLR_data$TPM)),"Gene"], 100))
      
      control_data <- rbind(control_data, 
                             PlantExpData[which(PlantExpData$Gene %in% candidateGenes),]) 
    }
  } else if (normalised == FALSE) {
    control_data <- PlantExpData[sample(nrow(PlantExpData[which(PlantExpData$TPM < max(NLR_data$TPM)),]), 1000),]
  }

  control_data$GeneSet <- rep("Control gene", times = nrow(control_data))
  
  
  # Merge the two datasets.
  PlantExpData <- control_data
  PlantExpData <- rbind(PlantExpData, NLR_data)
  
  # Add a column to PlantExpData with the expression level of each gene.
  expressionLevel <- c()
    
  for (row in 1:nrow(PlantExpData)) {
    if (PlantExpData[row, "TPM"] <= quantile(PlantExpData[PlantExpData$GeneSet=="R-gene",]$TPM, probs = .5)) {
      expressionLevel <- append(expressionLevel, "No Expression")
    }
    else if (PlantExpData[row, "TPM"] > quantile(PlantExpData[PlantExpData$GeneSet=="R-gene",]$TPM, probs = .5)) {
      expressionLevel <- append(expressionLevel, "Low Expression")
    }
  }
  
  PlantExpData <- cbind(PlantExpData, data.frame(ExpressionLevel = expressionLevel))
  
  # Sort control and R-genes into a hash based on expression level.
  sampleGenes <- hash()
  
  for (level in unique(PlantExpData$ExpressionLevel)) {
    for (set in c("R-gene", "Control gene")) {
      sampleGenes[[paste(set, level, sep = " ")]] <- PlantExpData[which(PlantExpData$Expression == level & PlantExpData$GeneSet == set),]
    }
  }
  return(sampleGenes)
}
