library(circlize)

brain.meth.dir <- c("../data/myDiff.brain.RData")
load(brain.meth.dir)


myDiff10p.hyper.brain <- getMethylDiff(myDiff.brain,
                                 difference=meth.cut,
                                 qvalue=qvalue.cut,
                                 type = "hyper")

myDiff10p.hypo.brain <- getMethylDiff(myDiff.brain,
                                       difference=meth.cut,
                                       qvalue=qvalue.cut,
                                       type = "hypo")

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

DMR_all$loc <- paste(DMR_all$chr, DMR_all$start, DMR_all$end, sep = "_")


annot <- t.brain %>% 
  filter(str_detect(external_gene_name, paste(rlv.genes.brain, collapse = "|"))) %>%
  tidyr::unite("loc", c("seqnames", "start", "end"), remove = FALSE) %>%
  dplyr::select(seqnames, start, end, mcols.meth.diff, external_gene_name, loc) %>%
  mutate(color = case_when(
    mcols.meth.diff >= 10 ~ "#660000",
    mcols.meth.diff <= -10 ~ "#000166"
  ),
  pointcolor = "green"
  )

pointcolor <- annot %>% dplyr::select(pointcolor, loc)

DMR_all <- left_join(DMR_all, pointcolor, by = "loc") %>%
  mutate(color = if_else(is.na(pointcolor), color, pointcolor))


tiff("./figures/circos_brain.tiff", 
     width = 2000, 
     height = 2000, 
     units = "px", 
     pointsize = 12,
     compression = "none",
     bg = "white",
     res = 300
) 


circos.clear()

circos.par(gap.after = c(rep(2, 9), 7), start.degree = 90)

circos.initializeWithIdeogram(plotType = c("labels", "axis"),
                              species = "xenTro10"
)
circos.par("track.height" = 0.3)


circos.genomicTrack(DMR_all, 
                    stack = FALSE, 
                    ylim = c(-100, 100),
                    panel.fun = function(region, value, ...) {
                      cex = 0.5
                      i = DMR_all[rownames(value), 5]
                      for(h in seq(-100, 100, by = 25)) {
                        circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
                      }
                      circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
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
circos.genomicLabels(annot, labels.column = 5, cex = 0.8, col= annot$color, line_lwd = 1, line_col="grey10", 
                                           side="inside", connection_height=0.05, labels_height=0.04)
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



