library(ggplot2)
library(data.table)
library(grid)
library(readr)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(rstudioapi)

allResultsProportions$GeneSet <- paste(str_match(allResultsProportions$GeneSet, "^([A-Za-z]+.gene).*$")[,-1], " \n", 
                                       str_match(allResultsProportions$GeneSet, "^[A-Za-z]+.gene(.*)$")[,-1], sep = "")


allResultsProportions <- allResultsProportions[order(factor(allResultsProportions$GeneSet, 
                                                              levels = c("Control gene \n No Expression", "R-gene \n No Expression",
                                                                         "Control gene \n Low Expression", "R-gene \n Low Expression"))),]
  
  
my_comparisons <- list(c("Control gene \n No Expression", "R-gene \n No Expression"), 
                       c("Control gene \n Low Expression", "R-gene \n Low Expression"),
                       c("R-gene \n No Expression", "R-gene \n Low Expression"))

for (size in unique(allResultsProportions$Size)) {
  if (size != "size") {
    df <- allResultsProportions[allResultsProportions$Size==size,]
    
    for (mod in unique(df$Mod.TF)) {
      df1 <- df[df$Mod.TF==mod,]
      
      mean_df <- data.frame()
      
      for (region in unique(df1$Region)) {
        df2 <- df1[df1$Region==region,]
        
        for (set in unique(df2$GeneSet)) {
          df3 <- df2[df2$GeneSet==set,]
          
          mean_df <- rbind(mean_df, data.frame(Size = size,
                                               Region = region,
                                               Mod.TF = mod,
                                               GeneSet = set,
                                               Enrichment.mean = mean(df3$Proportion),
                                               Enrichment.variance = sd(df3$Proportion),
                                               axisGroup = df3$axisGroup[1]))
        }
      }
      stat.test <- df1 %>% group_by(axisGroup) %>% 
        t_test(Proportion ~ GeneSet, comparisons = my_comparisons) %>% 
        mutate(y.position = rep(c(0.98, 0.98, 1.06), times = 10))
      
      plot <- ggbarplot(mean_df, x = "GeneSet", y="Enrichment.mean", ylab = "Average enrichment",
                        color = "black", fill = "GeneSet",
                        palette = c("azure3", "cadetblue", "bisque2", "lightsalmon2"), 
                        title = paste(mod, "-", length(unique(df2[which(grepl("R-gene", df2$GeneSet)),"Gene"])))) + theme_bw() +
        stat_pvalue_manual(
          stat.test, 
          label = "p.adj.signif", size = 4,
          tip.length = 0.01, hide.ns = FALSE) +

        coord_cartesian(ylim= c(0,1), clip = "off") +
        
        font("ylab", size = 14) +
        font("legend.title", size = 12) +
        font("legend.text", size = 10) +
        font("caption", size = 12) 
      
      plot <- facet(plot, facet.by = "axisGroup", nrow = 1, panel.labs.font = list(size = 10),
                    panel.labs = list(axisGroup = c("Intergenic","Promotor \n(1kb)","Promotor \n(500bp)", "20%",            
                                                    "40%","60%","80%","100%","Downstream \n(200bp)","Intergenic")))
      
      plot <- ggpar(plot, font.xtickslab = FALSE, ticks = FALSE, legend = "none", xlab = FALSE, legend.title = "",
                    font.ytickslab = 8)
      
      ggsave(paste("Graphs\\Enrichment\\PlantExp data\\", size, "-genes ", mod, ".png", sep = ""), plot = plot, width = 10, height = 4)  
    }
  } else next
}
