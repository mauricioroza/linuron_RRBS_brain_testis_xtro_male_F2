library(tidyverse)
library(svglite)
library(methylKit)
library(ggpubr)

gsea.brain <- readRDS("data/go.or.brain.rds")
gsea.brain.plot <- gsea.brain + 
  labs(title = "Gene Ontology Over-representation Analysis", subtitle = "Brain") +
  scale_y_discrete(labels = function(label) str_wrap(label, width = 25)) +
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,0.05), labels=c(0,0.025,0.05),breaks=c(0,0.025,0.05))



svglite("./figures/gsea_brain.svg",
        width = 15,
        height = 10,
        scaling = 2)

gsea.brain.plot

dev.off()

library(readxl)

df.kk.brain <- read_excel("supplementary_material_tables/S4_KEGG_enrichment.xlsx") %>%
  mutate_if(is.numeric, round, digits = 3)

kk.barplot <- ggplot(df.kk.brain, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,0.05), labels=c(0,0.025,0.05),breaks=c(0,0.025,0.05)) +
  labs(
    title = "KEGG Pathways Over-representation",
    subtitle = "Brain",
    x = "Count",
    y = NULL,
    fill = "p.adjust"
  ) + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, hjust = 0, colour = "black"),  # Align labels to the left
    axis.title = element_text(size = 12, colour = "black")
  ) +
  scale_y_discrete(labels = function(label) str_wrap(label, width = 25)) +  # Adjust the width as needed
  coord_fixed(ratio = (max(df.kk.brain$Count)/nrow(df.kk.brain)))

svglite("./figures/kk_brain.svg",
        width = 15,
        height = 10,
        scaling = 2)

kk.barplot

dev.off()

#horizontal legend

svglite("./figures/horizontal_legend.svg",
        width = 15,
        height = 10,
        scaling = 2)

kk.barplot + theme(legend.position = "bottom")

dev.off()
###############################################

gsea.testis <- readRDS("data/go.or.testis.rds")

gsea.testis.plot <- gsea.testis + labs(title = "Gene Ontology Over-representation Analysis", subtitle = "Testis")+
  scale_fill_gradient(low = "blue", high = "red", limits=c(0,0.05), labels=c(0,0.025,0.05),breaks=c(0,0.025,0.05))

svglite("./figures/gsea_testis.svg",
        width = 15,
        height = 10,
        scaling = 2)

gsea.testis.plot

dev.off()
