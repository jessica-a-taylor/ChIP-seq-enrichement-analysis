# Get expression data.
library(readxl)
library(readr)

# Import coordinates of the genomic regions of interest.
genomicData <- as.data.frame(read_csv("Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Import list of R-genes.
NLR_genes <- as.data.frame(read_xlsx("Arabidopsis NLRs.xlsx", sheet = 1))
NLR_genes <- genomicData[which(genomicData$Gene %in% NLR_genes$Gene),]

# Create a list of files in the folder for a particular tissue type.
filenames <- list.files(path = "PlantExp data\\Control\\" ,pattern="*.tsv")
expressionData <- data.frame()

# Merge the data from all files.
for (file in filenames) {
  expressionData <- rbind(expressionData, as.data.frame(read_tsv(paste("PlantExp data\\Control\\", file, sep = ""), show_col_types = FALSE)))
}

PlantExpData <- c()

# For each gene, calculate the mean expression across experiments.
for (row in 1:nrow(genomicData)) {
  PlantExpData <- rbind(PlantExpData, data.frame(Gene = genomicData[row,"Gene"],
                                                 seqnames = genomicData[row,"seqnames"],
                                                 start = genomicData[row,"start"],
                                                 end = genomicData[row,"end"],
                                                 width = genomicData[row,"width"],
                                                 strand = genomicData[row,"strand"],
                                                 ranges = genomicData[row,"ranges"],
                                                 TPM = mean(expressionData[expressionData$geneId == genomicData[row,"Gene"],"TPM"])))
}

write.csv(geneExpression, paste("PlantExp data\\geneExpression.csv", sep = "")) 