library(TxDb.Athaliana.BioMart.plantsmart28)
library(readr)
library(stringr)
library(hash)
library(readxl)


# Import all Arabidopsis genes.
Atgenes <- as.data.frame(transcriptsBy(TxDb.Athaliana.BioMart.plantsmart28, by="gene"))
colnames(Atgenes)[2] <- "Gene"

# Remove mitochondrial and chloroplast genes.
Atgenes$seqnames <- as.character(Atgenes$seqnames)
Atgenes <- Atgenes[-c(which(grepl("Mt", Atgenes$seqnames) | grepl("Pt", Atgenes$seqnames))),]

# Remove all non-coding genes.
geneTypes <- as.data.frame(read.table("Arabidopsis gene types.txt"))
geneTypes <- geneTypes[c(which(geneTypes$V2=="protein_coding")),]

Atgenes <- Atgenes[c(which(Atgenes$Gene %in% geneTypes$V1)),]

# Remove duplicate genes (different versions).
Atgenes <- Atgenes[c(which(Atgenes$tx_name == str_match(Atgenes$tx_name, "^([0-9a-zA-Z]+)([.])(1)$")[,1])),]

# Add 'ranges' column.
source("Functions\\Get range - merge gene coordinates.R")
Atgenes$ranges <- mergeCoordinates(Atgenes) 


# Assign each gene to euchromatic and heterochromatic regions.
source("Functions\\Overlaps functions.R")

pericentromericgeneRegions <- data.frame(Chromosome = c(1:5),
                                         Start = c("11500000", "1100000", "10300000", "1500000", "9000000"),
                                         End = c("17700000", "7200000", "17300000", "6300000", "16000000"))

Atgenes$Chromosomal.Region <- rep("", times = nrow(Atgenes))

for (row in 1:nrow(pericentromericgeneRegions)) {
  Atgenes[c(which(Atgenes$seqnames==row & 
                    betweenFunction(Atgenes$start, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==FALSE &
                    betweenFunction(Atgenes$end, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==FALSE)),"Chromosomal.Region"] <- "euchromatic"
  
  Atgenes[c(which(Atgenes$seqnames==row & 
                    betweenFunction(Atgenes$start, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==TRUE &
                    betweenFunction(Atgenes$end, as.numeric(pericentromericgeneRegions[row, "Start"]), as.numeric(pericentromericgeneRegions[row, "End"]))==TRUE)),"Chromosomal.Region"] <- "heterochromatic"
}

rm(pericentromericgeneRegions)

write.csv(Atgenes, file = "Protein coding genes.csv")