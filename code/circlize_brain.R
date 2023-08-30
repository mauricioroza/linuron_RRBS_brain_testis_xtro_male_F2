library(circlize)
library(tidyverse)
library(svglite)
library(methylKit)
  
brain.meth.dir <- c("./data/myDiff.brain.RData")
load(brain.meth.dir)

qvalue.cut <- 0.05
meth.cut <- 10

myDiff10p.brain <- getMethylDiff(myDiff.brain,
                                   difference=meth.cut,
                                   qvalue=qvalue.cut)
myDiff10p.brain <- myDiff10p.brain[order(myDiff10p.brain$qvalue),]
  
myDiff10p.hyper.brain <- getMethylDiff(myDiff.brain,
                                         difference=meth.cut,
                                         qvalue=qvalue.cut,
                                         type = "hyper")
  
myDiff10p.hypo.brain <- getMethylDiff(myDiff.brain,
                                        difference=meth.cut,
                                        qvalue=qvalue.cut,
                                        type = "hypo")
  
t.brain <- readRDS("./data/t.brain.rds")
  


DMR_hyper <- data.frame("chr" = myDiff10p.hyper.brain$chr,
                        "start" = myDiff10p.hyper.brain$start,
                        "end" = myDiff10p.hyper.brain$end,
                        "meth.diff" = myDiff10p.hyper.brain$meth.diff)

DMR_hypo <- data.frame("chr" = myDiff10p.hypo.brain$chr,
                       "start" = myDiff10p.hypo.brain$start,
                       "end" = myDiff10p.hypo.brain$end,
                       "meth.diff" = myDiff10p.hypo.brain$meth.diff)

bed_list <- list(DMR_hyper, DMR_hypo)


####All DMRs

all.DMRs <- data.frame(myDiff.brain) %>% filter(qvalue <= 0.05)

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



rlv.genes.brain <- c("grik2", "slc17a7", "grm4", "grm5", "grin1", "grm5", "grm1", "grm8", "grin2b", #glutamate signaling
                     "gabbr1", "gabbr2", "gabrb3", #GABA signaling
                     "dnmt3a", "hdac8", "mbd2", #methylation
                     "gh1", "irs2", "igfbp4", "igfbp5", #somatotropic
                     "trhr3", "trhde", "dio1", "\\btg\\b", #thyrotropic 
                     "esrrg",
                     "kiss2", "prkag2", "prkar1b", #GnRH signaling
                     "hsd17b12") #steroidogenesis

one <- c("dnmt3a", "hdac8", "mbd2") #methylation
two <- c("kiss2", "prkag2", "prkar1b") #GnRH signaling
three <- c("esrrg", "hsd17b12") #steroidogenesis
four <- c("gh1", "irs2", "igfbp4", "igfbp5") #somatotropic
five <- c("trhr3", "trhde", "dio1", "\\btg\\b") #thyrotropic
six <- c("gabbr1", "gabbr2", "gabrb3") #GABA signaling
seven <- c("grik2", "slc17a7", "grm4", "grm5", "grin1", "grm5", "grm1", "grm8", "grin2b") #glutamate signaling



hyper_color_rlv <- "#660000"
hypo_color_rlv <- "#000166"

annot <- t.brain %>%
  filter(str_detect(external_gene_name, paste(rlv.genes.brain, collapse = "|"))) %>%
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
                                                     "³\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(four, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(four, collapse = "|"), 
                                                     "⁴\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(five, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(five, collapse = "|"), 
                                                     "⁵\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(six, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(six, collapse = "|"), 
                                                     "⁶\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(seven, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(seven, collapse = "|"), 
                                                     "⁷\\0"),
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


circos.brain <- function() {
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
circos.genomicLabels(annot, labels.column = 5, cex = 0.8, col= annot$pointcolor, line_lwd = 1, line_col="grey10", 
                                           side="inside", connection_height=0.05, labels_height=0.04)

}
tiff("./figures/circos_brain.tiff", 
     width = 2000, 
     height = 2000, 
     units = "px", 
     pointsize = 12,
     compression = "none",
     bg = "white",
     res = 300
) 

circos.brain()

dev.off()

svglite("./figures/circos_brain.svg",
        width = 18,
        height = 18,
        scaling = 2.2)

circos.brain()

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



