RNA_seqAnalysis <- function(dataToUse, exLevel) {
  
  treatments <- c("Mock", "untreated", "control", "Not treated", "mock", "unstressed control", "none",
                  "normal condition","healthy", "NA")
  
  tissue <- c("leaves", "seedlings")
  
  expressionData <- hash()
  
  for (test in names(dataToUse)) {
    print(test)
    
    expressionData[[test]] <- as.data.frame(read_csv(paste("RNA-seq data\\", test, "_genes.csv", sep = ""), skip = 1))
    
    expressionData[[test]] <- expressionData[[test]][,-c(1:2, ncol(expressionData[[test]]):(ncol(expressionData[[test]])-3))]
    
    expressionData[[test]] <- expressionData[[test]][expressionData[[test]]$Ecotype=="Col-0",]
    expressionData[[test]] <- expressionData[[test]][expressionData[[test]]$Genotype=="wild type",]
    
    expressionData[[test]] <- expressionData[[test]][c(which(expressionData[[test]]$Treatment %in% treatments)),]
    expressionData[[test]] <- expressionData[[test]][c(which(expressionData[[test]]$Tissue %in% tissue)),]
    
    meanExpressionData <- data.frame()
    
   for (col in 1:ncol(expressionData[[test]][,c(which(!is.na(str_match(colnames(expressionData[[test]]), "^AT[1-9]+G[0-9]+$"))))])) {
      meanExpressionData <- rbind(meanExpressionData, data.frame(Gene = colnames(expressionData[[test]])[col],
                                                                         Expression = mean(expressionData[[test]][,col])))
      }
      
    meanExpressionData <- meanExpressionData[order(meanExpressionData$Gene),]
      
    expressionLevel <- c()
    
    for (gene in meanExpressionData$Gene) {
      if (0 <= meanExpressionData[meanExpressionData$Gene==gene, "Expression"] & meanExpressionData[meanExpressionData$Gene==gene, "Expression"] <= 1) {
        expressionLevel <- append(expressionLevel, "No Expression")
      }
      else if (1 < meanExpressionData[meanExpressionData$Gene==gene, "Expression"] & meanExpressionData[meanExpressionData$Gene==gene, "Expression"] <= 50) {
        expressionLevel <- append(expressionLevel, "Low Expression")
      }
      else if (50 < meanExpressionData[meanExpressionData$Gene==gene, "Expression"] & meanExpressionData[meanExpressionData$Gene==gene, "Expression"] <= 100) {
        expressionLevel <- append(expressionLevel, "Intermediate Expression")
      }
      #else if (50 < meanExpressionData[meanExpressionData$Gene==gene, "Expression"] & meanExpressionData[meanExpressionData$Gene==gene, "Expression"]<= 100) {
       # expressionLevel <- append(expressionLevel, "High Expression")
      #}
      else if (meanExpressionData[meanExpressionData$Gene==gene, "Expression"] > 100) {
        expressionLevel <- append(expressionLevel, "High Expression")
      }
    }
    
    meanExpressionData <- cbind(meanExpressionData, data.frame(ExpressionLevel = expressionLevel))
    # Sort genes to hashes based on tissue and expression level.
   for (level in exLevel) {
      df1 <- meanExpressionData[meanExpressionData$ExpressionLevel==level,]
      df1 <- df1[c(which(df1$Gene %in% dataToUse[[test]]$Gene)),]
      
      dataToUse[[paste(test, "_", level, sep = "")]] <- dataToUse[[test]][c(which(dataToUse[[test]]$Gene %in% df1$Gene)),]
      
      dataToUse[[paste(test, "_", level, sep = "")]] <- cbind(dataToUse[[paste(test, "_", level, sep = "")]], df1[,c(2:3)])
    }   
    write.csv(meanExpressionData, file = "RNA-seq expression data.csv")
  }
  return(dataToUse)
}
