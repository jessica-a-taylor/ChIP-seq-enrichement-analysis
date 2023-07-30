for (normalised in c(TRUE, FALSE)) {
  if (normalised == FALSE) {
    allResultsFrequencies <- data.frame(read_csv(paste(analysis, "\\Non-normalised\\allResultsFrequencies.csv", sep = "")))
  } else if (normalised == TRUE) {
    allResultsFrequencies <- data.frame(read_csv(paste(analysis, "\\Normalised\\allResultsFrequencies.csv", sep = "")))
  }
  
  # Replace comma in 'Comparisons' column with \n.
  allResultsFrequencies$Comparison <- gsub(",", "\n", allResultsFrequencies$Comparison)
  
  for (mod in unique(allResultsFrequencies$Mod.TF)) {
    df <- allResultsFrequencies[allResultsFrequencies$Mod.TF==mod,]
    
    plot <- ggbarplot(df, x = "Comparison", y="Enrichment.variance", ylab = "Enrichment variance",
                      color = "black", fill = "Comparison", 
                      palette = c("azure3", "cadetblue", "bisque2", "darksalmon"), title = mod) + 
      theme_bw() +
      
      font("title", size = 16) +
      font("ylab", size = 14) +
      font("legend.title", size = 14) +
      font("legend.text", size = 12) +
      font("caption", size = 12) 
    
    plot <- facet(plot, facet.by = "axisGroup", nrow = 1, panel.labs.font = list(size = 10),
                  panel.labs = list(axisGroup = c("Intergenic","Promotor \n(1kb)","Promotor \n(500bp)", "20%",            
                                                  "40%","60%","80%","100%","Downstream \n(200bp)","Intergenic")))
    
    plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "bottom", xlab = FALSE, legend.title = "",
                  font.ytickslab = 8)
    
    if (normalised == FALSE) {
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\Non-normalised\\", mod, "_variance.png", sep = ""), plot = plot, width = 10, height = 4)  
    } else if (normalised == TRUE) {
      ggsave(paste("Graphs\\Enrichment\\", analysis, "\\Normalised\\", mod, "_variance.png", sep = ""), plot = plot, width = 10, height = 4)  
    }
  }
}