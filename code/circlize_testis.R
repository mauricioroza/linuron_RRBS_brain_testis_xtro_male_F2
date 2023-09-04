library(circlize)
library(tidyverse)
library(svglite)
library(methylKit)

testis.meth.dir <- c("./data/myDiff.testis.RData")
load(testis.meth.dir)

qvalue.cut <- 0.05
meth.cut <- 10

myDiff10p.testis <- getMethylDiff(myDiff.testis,
                                 difference=meth.cut,
                                 qvalue=qvalue.cut)
myDiff10p.testis <- myDiff10p.testis[order(myDiff10p.testis$qvalue),]

myDiff10p.hyper.testis <- getMethylDiff(myDiff.testis,
                                       difference=meth.cut,
                                       qvalue=qvalue.cut,
                                       type = "hyper")

myDiff10p.hypo.testis <- getMethylDiff(myDiff.testis,
                                      difference=meth.cut,
                                      qvalue=qvalue.cut,
                                      type = "hypo")

t.testis <- readRDS("./data/t.testis.rds")



DMR_hyper <- data.frame("chr" = myDiff10p.hyper.testis$chr,
                        "start" = myDiff10p.hyper.testis$start,
                        "end" = myDiff10p.hyper.testis$end,
                        "meth.diff" = myDiff10p.hyper.testis$meth.diff)

DMR_hypo <- data.frame("chr" = myDiff10p.hypo.testis$chr,
                       "start" = myDiff10p.hypo.testis$start,
                       "end" = myDiff10p.hypo.testis$end,
                       "meth.diff" = myDiff10p.hypo.testis$meth.diff)

bed_list <- list(DMR_hyper, DMR_hypo)


####All DMRs

all.DMRs <- data.frame(myDiff.testis) %>% filter(qvalue <= 0.05)

DMR_all <- data.frame("chr" = all.DMRs$chr,
                      "start" = all.DMRs$start,
                      "end" = all.DMRs$end,
                      "meth.diff" = all.DMRs$meth.diff)

DMR_all <- DMR_all %>%
  mutate(color = case_when(
    meth.diff >= 10 ~ "#E41A1C",
    meth.diff <= -10 ~ "#377EB8",
    meth.diff < 10 & meth.diff > -10 ~ "#D3D3D3"
  )) 



rlv.genes.testis <- c("ep300", "elp3", "kat5", "kat14", #histone acetyltransferase
                      "lhcgr", "hsd17b12", "esrrg",
                      "piwil1", "mael", "spo11", "\\bddx4\\b", #spermatogenesis, gonadal development
                      "dnmt3a"
)

one <- c("dnmt3a", "ep300", "elp3", "kat5", "kat14") #methylation, epigenetic mechanisms
two <- c("lhcgr") #GnRH signalling
three <- c("hsd17b12", "esrrg", "lhcgr")  #steroidogenesis
eight <- c("piwil1", "mael", "spo11", "\\bddx4\\b") #spermatogenesis, gonadal development
nine <- c("ep300", "elp3", "kat5", "kat14") #histone acetyltransferase


hyper_color_rlv <- "#660000"
hypo_color_rlv <- "#000166"

annot <- t.testis %>%
  filter(str_detect(external_gene_name, paste(rlv.genes.testis, collapse = "|"))) %>%
  tidyr::unite("loc", c("seqnames", "start", "end"), remove = FALSE) %>%
  dplyr::select(seqnames, start, end, mcols.meth.diff, external_gene_name, loc) %>%
  mutate(pointcolor = case_when(
    mcols.meth.diff >= 10 ~ hyper_color_rlv,
    mcols.meth.diff <= -10 ~ hypo_color_rlv
  )) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(one, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(one, collapse = "|"), 
                                                     "¹\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(two, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(two, collapse = "|"), 
                                                     "²\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(three, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(three, collapse = "|"), 
                                                     "⁸\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(eight, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(eight, collapse = "|"), 
                                                     "⁸\\0"),
                                     external_gene_name))

pointcolor <- annot %>% dplyr::select(pointcolor, loc)

DMR_all2 <- DMR_all
DMR_all2$loc <- paste(DMR_all$chr, DMR_all$start, DMR_all$end, sep = "_")

DMR_all2 <- left_join(DMR_all2, pointcolor, by = "loc") %>%
  mutate(color = if_else(is.na(pointcolor), color, pointcolor),
         cex = case_when(
           color == hyper_color_rlv ~ 1,
           color == hypo_color_rlv ~ 1,
           TRUE ~ 0.5
         )) %>%
  dplyr::select(-loc, -pointcolor)

DMR_all2 <- DMR_all2[order(DMR_all2$cex),]


circos.testis <- function() {
  circos.clear()
  
  circos.par(gap.after = c(rep(2, 9), 7), start.degree = 90)
  
  circos.initializeWithIdeogram(plotType = c("labels", "axis"),
                                species = "xenTro10"
  )
  circos.par("track.height" = 0.3)
  
  
  circos.genomicTrack(DMR_all2[,1:4], 
                      stack = FALSE, 
                      ylim = c(-100, 100),
                      panel.fun = function(region, value, ...) {
                        cex = DMR_all2[rownames(value), "cex"]
                        i = DMR_all2[rownames(value), "color"]
                        for(h in seq(-100, 100, by = 25)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
                        }
                        circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, bg = i, ...)
                        circos.lines(CELL_META$cell.xlim, c(-10, -10), col = "#7D7D7D", lwd = 1, lty = "dashed")
                        circos.lines(CELL_META$cell.xlim, c(10, 10), col = "#7D7D7D", lwd = 1, lty = "dashed")
                        circos.lines(CELL_META$cell.xlim, c(0, 0), col = "#888888", lwd = 1, lty = "dotted")
                        circos.yaxis(side = "left", 
                                     at = c(-100, -50, 0, 50, 100),
                                     labels = paste0(c(-100, -50, 0, 50, 100), "%"),
                                     sector.index = get.all.sector.index()[1], 
                                     labels.cex = 0.3
                        )
                      })  
  circos.genomicLabels(annot, labels.column = 5, cex = 1, col= annot$pointcolor, line_lwd = 1, line_col="grey10", 
                       side="inside", connection_height=0.025, labels_height=0.05, padding = 0.8)
  
}
tiff("./figures/circos_testis.tiff", 
     width = 2000, 
     height = 2000, 
     units = "px", 
     pointsize = 12,
     compression = "none",
     bg = "white",
     res = 300
) 

circos.testis()

dev.off()

svglite("./figures/circos_testis.svg",
        width = 16,
        height = 16,
        scaling = 2)

circos.testis()

dev.off()



# circos.clear()
# 
# circos.par(gap.after = c(rep(2, 9), 7), start.degree = 90)
# 
# circos.initializeWithIdeogram(plotType = c("labels", "axis"),
#                               species = "xenTro10"
# )
# circos.par("track.height" = 0.3)
# 
# circos.genomicTrack(DMR_all, 
#                     stack = FALSE, 
#                     ylim = c(-100, 100),
#                     panel.fun = function(region, value, ...) {
#                       cex = 0.5
#                       i = ifelse(value[[1]] > 10, "#E41A1C", ifelse(value[[1]] > -10, "#D3D3D3", "#377EB8"))
#                       for(h in seq(-100, 100, by = 25)) {
#                         circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
#                       }
#                       circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
#                       circos.lines(CELL_META$cell.xlim, c(-10, -10), col = "#7D7D7D", lwd = 1, lty = "dashed")
#                       circos.lines(CELL_META$cell.xlim, c(10, 10), col = "#7D7D7D", lwd = 1, lty = "dashed")
#                       circos.lines(CELL_META$cell.xlim, c(0, 0), col = "#888888", lwd = 1, lty = "dotted")
#                       circos.yaxis(side = "left", 
#                                    at = c(-100, -50, 0, 50, 100),
#                                    labels = paste0(c(-100, -50, 0, 50, 100), "%"),
#                                    sector.index = get.all.sector.index()[1], 
#                                    labels.cex = 0.3
#                       )
#                     })  
# 
# circos.par("track.height" = 0.075)
# 
# circos.genomicDensity(bed_list, col = c("#FF000080", "#0000FF80"))
# 
# circos.genomicLabels(annot, labels.column = 5, cex = 0.8, col= annot$color, line_lwd = 1, line_col="grey10", 
#                      side="inside", connection_height=0.05, labels_height=0.04)
# 
# dev.off()
# 
# 
