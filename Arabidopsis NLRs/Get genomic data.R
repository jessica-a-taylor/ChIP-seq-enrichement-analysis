# Import genome.
genomeData <- as.data.frame(read.csv("Arabidopsis NLRs/Data/Arabidopsis_thaliana.TAIR10.cds.all.csv"))
genomeData <- genomeData[order(genomeData$GeneID),]

# Remove mitochondrial and chloroplast genes.
genomeData <- genomeData[-c(which(grepl("Mt", genomeData$seqnames) | grepl("Pt", genomeData$seqnames))),]

# Remove duplicate genes.
genomeData <- genomeData[which(genomeData$GeneID == str_match(genomeData$GeneID, "^([0-9a-zA-Z]+)([.])(1)$")[,1]),]

# Add 'width' and 'ranges' columns.
source("Functions/Gene width&range.R")
genomeData$width <- getWidth(genomeData)
genomeData$ranges <- getRange(genomeData)

# Change the '+/-' strand to 'positive/negative'.
genomeData[which(genomeData$strand=="1"), "strand"] <- "positive"
genomeData[which(genomeData$strand=="-1"), "strand"] <- "negative"

# Assign each gene to euchromatic and heterochromatic regions.
# Create function that determines whether value a is between values b and c.
betweenFunction <- function(a,b,c) {
  return(b<a & a<c)
}

pericentromericgeneRegions <- data.frame(start = c("11500000", "1100000", "10300000", "1500000", "9000000"),
                                         end = c("17700000", "7200000", "17300000", "6300000", "16000000"))

genomeData$Chromosomal.Region <- rep("", times = nrow(genomeData))

for (row in 1:nrow(pericentromericgeneRegions)) {
  genomeData[c(which(genomeData$seqnames==row & 
                    betweenFunction(genomeData$start, as.numeric(pericentromericgeneRegions[row, "start"]), as.numeric(pericentromericgeneRegions[row, "end"]))==FALSE &
                    betweenFunction(genomeData$end, as.numeric(pericentromericgeneRegions[row, "start"]), as.numeric(pericentromericgeneRegions[row, "end"]))==FALSE)),"Chromosomal.Region"] <- "euchromatic"
  
  genomeData[c(which(genomeData$seqnames==row & 
                    betweenFunction(genomeData$start, as.numeric(pericentromericgeneRegions[row, "start"]), as.numeric(pericentromericgeneRegions[row, "end"]))==TRUE &
                    betweenFunction(genomeData$end, as.numeric(pericentromericgeneRegions[row, "start"]), as.numeric(pericentromericgeneRegions[row, "end"]))==TRUE)),"Chromosomal.Region"] <- "heterochromatic"
}

rm(pericentromericgeneRegions)

write.csv(genomeData, file = "Arabidopsis NLRs/Data/Protein coding genes.csv")