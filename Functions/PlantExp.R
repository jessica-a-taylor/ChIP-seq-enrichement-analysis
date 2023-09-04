# Function to get filtered expression data for each set of sample genes in each tissue. 
PlantExp <- function(normalised) {
  
  # Import expression data.
  PlantExpData <- as.data.frame(read_csv("PlantExp data\\geneExpression.csv"))
  PlantExpData <- PlantExpData[-c(which(is.na(PlantExpData$TPM))),]
  
  # Import list of R-genes.
  NLR_data <- as.data.frame(read_xlsx("Arabidopsis NLRs.xlsx", sheet = 1))
  NLR_data <- PlantExpData[c(which(PlantExpData$Gene %in% NLR_data$Gene)),]
  NLR_data$GeneSet <- rep("R-gene", times = nrow(NLR_data))
  
  # Sample 1000 random genes with a gene length distribution equal to that of the R-genes.
  if (normalised == TRUE) {
    euchromaticGenes <- c()
    heterochromaticGenes <- c()
    
    for (n in seq(from = .1, to = 1, by = .1)) {
      euchromaticGenes <- append(euchromaticGenes, sample(PlantExpData[which(PlantExpData$width > quantile(NLR_data$width, probs = n-.1) &
                                                                          PlantExpData$width <= quantile(NLR_data$width, probs = n) &
                                                                          !(PlantExpData$Gene %in% euchromaticGenes) &
                                                                          !(PlantExpData$Gene %in% NLR_data$Gene) &
                                                                           PlantExpData$TPM < max(NLR_data$TPM)),"Gene"], 
                                                      100*(nrow(NLR_data[which(NLR_data$width > quantile(NLR_data$width, probs = n-.1) &
                                                                               NLR_data$width <= quantile(NLR_data$width, probs = n) &
                                                                               NLR_data$Chromosomal.Region=="euchromatic"),])/nrow(NLR_data[which(NLR_data$width > quantile(NLR_data$width, probs = n-.1) &
                                                                                                                                                   NLR_data$width <= quantile(NLR_data$width, probs = n)),]))))
      
      heterochromaticGenes <- append(heterochromaticGenes, sample(PlantExpData[which(PlantExpData$width > quantile(NLR_data$width, probs = n-.1) &
                                                                             PlantExpData$width <= quantile(NLR_data$width, probs = n) &
                                                                             !(PlantExpData$Gene %in% heterochromaticGenes) &
                                                                             !(PlantExpData$Gene %in% NLR_data$Gene) &
                                                                             PlantExpData$TPM < max(NLR_data$TPM)),"Gene"], 
                                                        100*(nrow(NLR_data[which(NLR_data$width > quantile(NLR_data$width, probs = n-.1) &
                                                                                 NLR_data$width <= quantile(NLR_data$width, probs = n) &
                                                                                 NLR_data$Chromosomal.Region=="heterochromatic"),])/nrow(NLR_data[which(NLR_data$width > quantile(NLR_data$width, probs = n-.1) &
                                                                                                                                                         NLR_data$width <= quantile(NLR_data$width, probs = n)),]))))
      
    }
    control_data <- PlantExpData[which(PlantExpData$Gene %in% c(euchromaticGenes, heterochromaticGenes)),]
    
  } else if (normalised == FALSE) {
    control_data <- PlantExpData[sample(nrow(PlantExpData(which(PlantExpData$TPM < max(NLR_data$TPM)))), 1000),]
  }

  control_data$GeneSet <- rep("Control gene", times = nrow(control_data))
  
  # Merge the two datasets.
  PlantExpData <- control_data
  PlantExpData <- rbind(PlantExpData, NLR_data)
  PlantExpData <- PlantExpData[order(PlantExpData$Gene),]

  
  # Add a column to PlantExpData with the expression level of each gene.
  expressionLevel <- c()
  
  for (row in 1:nrow(PlantExpData)) {
    if (PlantExpData[row, "TPM"] <= 5) {
      expressionLevel <- append(expressionLevel, "No Expression")
    }
    else if (PlantExpData[row, "TPM"] > 5) {
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