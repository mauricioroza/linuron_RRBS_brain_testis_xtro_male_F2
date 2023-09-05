library(tidyverse)
library(svglite)
library(methylKit)
library(ggpubr)

gsea.brain <- readRDS("data/go.or.brain.rds")
gsea.brain[["data"]] <- gsea.brain[["data"]][c(1:3,5,6),]
gsea.brain.plot <- gsea.brain + 
  labs(title = "Gene Ontology\nOver-representation", subtitle = "Brain") +
  scale_y_discrete(labels = function(label) str_wrap(label, width = 25)) +
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,0.05), labels=c(0,0.025,0.05),breaks=c(0,0.025,0.05)) +
  theme(
    axis.text = element_text(size = 16),        # Adjust axis label text size
    axis.title = element_text(size = 16),       # Adjust axis title text size
    plot.title = element_text(size = 18),        # Adjust plot title text size
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  coord_fixed(ratio = (max(gsea.brain[["data"]]$Count)/nrow(gsea.brain[["data"]])))



svglite("./figures/gsea_brain.svg",
        width = 5,
        height = 5,
        scaling = 1)

gsea.brain.plot

dev.off()

library(readxl)

df.kk.brain <- read_excel("supplementary_material_tables/S4_KEGG_enrichment.xlsx") %>%
  mutate_if(is.numeric, round, digits = 3)

kk.barplot <- ggplot(df.kk.brain, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,0.05), labels=c(0,0.025,0.05),breaks=c(0,0.025,0.05)) +
  labs(
    title = "KEGG Pathways\nOver-representation",
    subtitle = "Brain",
    x = "Count",
    y = NULL,
    fill = "p.adjust"
  ) + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 15, hjust = 0, colour = "black"),  # Align labels to the left
    axis.title = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 18),        # Adjust plot title text size
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  scale_y_discrete(labels = function(label) str_wrap(label, width = 25)) +  # Adjust the width as needed
  coord_fixed(ratio = (max(df.kk.brain$Count)/nrow(df.kk.brain)))

svglite("./figures/kk_brain.svg",
        width = 5,
        height = 5,
        scaling = 1)

kk.barplot

dev.off()

#horizontal legend

svglite("./figures/horizontal_legend.svg",
        width = 5,
        height = 5,
        scaling = 1)

kk.barplot + theme(legend.position = "bottom")

dev.off()
###############################################

gsea.testis <- readRDS("data/go.or.testis.rds")

gsea.testis.plot <- gsea.testis + labs(title = "Gene Ontology\nOver-representation", subtitle = "Testis")+
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,0.05), labels=c(0,0.025,0.05),breaks=c(0,0.025,0.05)) +
  theme(
    axis.text = element_text(size = 16),        # Adjust axis label text size
    axis.title = element_text(size = 16),       # Adjust axis title text size
    plot.title = element_text(size = 18),        # Adjust plot title text size
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  coord_fixed(ratio = (max(gsea.testis[["data"]]$Count)/nrow(gsea.testis[["data"]])))

svglite("./figures/gsea_testis.svg",
        width = 5,
        height = 5,
        scaling = 1)

gsea.testis.plot

dev.off()
