# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(genomicData, NLR_genes) {

  # Create a list of files in the folder for a particular tissue type.
  filenames <- list.files(path = "PlantExp data\\Control\\" ,pattern="*.tsv")
  expressionData <- data.frame()
  
  # Merge the data from all files.
  for (file in filenames) {
    expressionData <- rbind(expressionData, as.data.frame(read_tsv(paste("PlantExp data\\Control\\", file, sep = ""), show_col_types = FALSE)))
  }
  
  # Sample 1000 random genes with a gene length distribution equal to that of the R-genes.
  control_genes <- data.frame()
  for (n in seq(from = .1, to = 1, by = .1)) {
    control_genes <- rbind(control_genes, 
                           genomicData[c(sample(nrow(genomicData[which(genomicData$width > quantile(NLR_genes$width, probs = n-.1) &
                                                                       genomicData$width <= quantile(NLR_genes$width, probs = n)),]), 100)),]) 
  }

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
    if (PlantExpData[row, "TPM"] <= quantile(PlantExpData[PlantExpData$GeneSet=="R-gene",]$TPM, probs = .5)) {
      expressionLevel <- append(expressionLevel, "No Expression")
    }
    #else if (quantile(PlantExpData[PlantExpData$GeneSet=="R-gene",]$TPM, probs = .5) < PlantExpData[row, "TPM"] & PlantExpData[row, "TPM"] <= quantile(PlantExpData[PlantExpData$GeneSet=="R-gene",]$TPM, probs = .8)) {
    #  expressionLevel <- append(expressionLevel, "Low Expression")
    #}
    else if (quantile(PlantExpData[PlantExpData$GeneSet=="R-gene",]$TPM, probs = .5) < PlantExpData[row, "TPM"]) {
      expressionLevel <- append(expressionLevel, "Low Expression")
    }
    #else if (PlantExpData[row, "TPM"] > 10) {
    #  expressionLevel <- append(expressionLevel, "High Expression")
    #}
  }
  
  PlantExpData <- cbind(PlantExpData, data.frame(ExpressionLevel = expressionLevel))
  
  # Add genomic data to the `PlantExpData` dataframe.
  infoNeeded <- data.frame()
  
  for (row in 1:nrow(PlantExpData)) {
    if (PlantExpData[row, "GeneSet"] == "Control gene") {
      infoNeeded <- rbind(infoNeeded, data.frame(control_genes[which(control_genes$Gene == PlantExpData[row, "Gene"]),][1,]))
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
      sampleGenes[[paste(set, level, sep = " ")]] <- PlantExpData[which(PlantExpData$Expression == level & PlantExpData$GeneSet == set),]
    }
  }
  return(sampleGenes)
}
