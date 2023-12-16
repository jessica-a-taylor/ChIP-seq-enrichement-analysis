library(rtracklayer)
library(glue)
library(stringr)

getChIP_seq_data <- function() {
  source("Functions/Gene width&range.R")
  
  ChIP_experiments <- as.data.frame(read.xlsx("Arabidopsis NLRs/Nextflow_backup/ChIP experiment SRA data.xlsx"))

  # Combine nextflow pipeline outputs into a single dataframe.
  nextflowOutput <- data.frame()
  
  for (file in list.files(path = "Arabidopsis NLRs/Nextflow_backup/", pattern = "Peaks.bed")) {
    
    # Rename file to include target modification/TF (if it has not been changed already).
    if (str_detect(file, "_SRR") == FALSE) {
      
      file.rename(paste("Arabidopsis NLRs/Nextflow_backup/", file, sep = ""), 
                  paste("Arabidopsis NLRs/Nextflow_backup/", ChIP_experiments[which(str_detect(ChIP_experiments$Sample.data, str_match(file,"^(SRR.*)_merged.*$")[,2])==TRUE), "Mod.TF"],
                        "_", file, sep = "")) 
      
      file <- paste(ChIP_experiments[which(str_detect(ChIP_experiments$Sample.data, str_match(file,"^(SRR.*)_merged.*$")[,2])==TRUE), "Mod.TF"],
                    "_", file, sep = "")
    }
    
    # Merge all broad and narrow peaks datasets.
    data <-  as.data.frame(import.bed(paste("Arabidopsis NLRs/Nextflow_backup/", file, sep = "")))
    
    # Add a column with the experiment code.
    data$experiment <- rep(str_match(file, "^.*_(SRR[0-9]+).*$")[,-1], times = nrow(data))
    
    # Add data to 'nextflowOutput'.
    nextflowOutput <- rbind(nextflowOutput, data)
  }
  
  nextflowOutput <- nextflowOutput[-which(nextflowOutput$seqnames == "ChrC" | nextflowOutput$seqnames == "ChrM"),]
  nextflowOutput$seqnames <- str_match(nextflowOutput$seqnames, "^Chr([0-9]+)$")[,2]
  
  # Add ranges column to 'nextflowOutput'.
  nextflowOutput$ranges <- getRange(nextflowOutput) 
  
  # Add column containing the chromatin modification/TF investigated.
  nextflowOutput$`Mod.TF` <- rep(NA, times = nrow(nextflowOutput))
  
  for (mod in unique(ChIP_experiments$`Modification/TF`)) {
    focusModification <- ChIP_experiments[ChIP_experiments$`Modification/TF`==mod,]
    
    # Create a list into which the ChIP-seq experiments for the focus modification/TF will be stored.
    focusExperiments <- c()
    
    for (row in 1:nrow(focusModification)) {
      focusExperiments <- append(focusExperiments, focusModification[row, "Sample.data"])
    }
    # Convert space-separated experiment list into a comma-separated list.
    focusExperiments <- glue_collapse(focusExperiments, ",")
    focusExperiments <- gsub(" ", ",", focusExperiments)
    focusExperiments <- c(strsplit(focusExperiments, ",")[[1]])
    
    # Replace 'NAs' in 'Mod.TF' column with the focus modification/TF
    nextflowOutput[which(nextflowOutput$experiment %in% focusExperiments),"Mod.TF"] <- mod
  }
  
  # Remove rows with any remaining 'NAs'.
  nextflowOutput <- nextflowOutput[c(which(!is.na(nextflowOutput$Mod.TF))),]
  
  return(nextflowOutput)
}