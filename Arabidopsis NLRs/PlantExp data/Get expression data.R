library(readr)

# Import coordinates of the genomic regions of interest.
genomicData <- as.data.frame(read.csv("Arabidopsis NLRs/Data/Protein coding genes.csv"))
genomicData <- genomicData[,-1]

# Create a list of files in the folder for a particular tissue type.
filenames <- list.files(path = "Arabidopsis NLRs/PlantExp data/Control conditions/" ,pattern="*.tsv")
expressionData <- data.frame()

# Merge the data from all files.
for (file in filenames) {
  expressionData <- rbind(expressionData, 
                          as.data.frame(read_tsv(paste("Arabidopsis NLRs/PlantExp data/Control conditions/", file, sep = ""), 
                                                 show_col_types = FALSE)))
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
                                                 Chromosomal.Region = genomicData[row, "Chromosomal.Region"],
                                                 TPM = mean(expressionData[expressionData$geneId == genomicData[row,"Gene"],"TPM"])))
}

write.csv(PlantExpData, paste("Arabidopsis NLRs/PlantExp data/geneExpression.csv", sep = "")) 