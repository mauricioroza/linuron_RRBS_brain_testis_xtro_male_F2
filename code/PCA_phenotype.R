library(tidyverse)
library(readxl)
library(ChIPpeakAnno)


testis.meth.dir <- c("./data/myDiff.testis.RData")
load(testis.meth.dir)



phenotype <- read_excel("data/phenotype_data.xlsx")

brain.unite <- readRDS("./data/brain_united.rds")

pm.brain <- percMethylation(brain.unite) %>% data.frame

rownames(pm.brain)

brain.meth.dir <- c("./data/myDiff.brain.RData")
load(brain.meth.dir)


meth.cut <- 10
qvalue.cut <- 0.05

diff.meth.brain <- getMethylDiff(myDiff.brain,
                                 difference=meth.cut,
                                 qvalue=qvalue.cut) %>% getData()

pm.brain.dm <- pm.brain[rownames(pm.brain) %in% rownames(diff.meth.brain),]

pm.brain.dm.merged <- cbind(pm.brain.dm, diff.meth.brain)

t.brain <- read_excel("./supplementary_material_tables/S1_DMR_and_overlapping_genes.xlsx", sheet = "Brain")


meth_pos.pca.brain <- GRanges(
  seqnames = pm.brain.dm.merged$chr, 
  ranges = IRanges(start = pm.brain.dm.merged$start, end = pm.brain.dm.merged$end),
  mcols = pm.brain.dm.merged[, c(1:16,21:23)]
)

martfile <- c("./data/mart_xtropicalis_gene_ensembl.RData")

load(martfile)

xtro <- getAnnotation(mart, featureType = "TSS")

pca.brain <- annotatePeakInBatch(meth_pos.pca.brain, AnnotationData = xtro, FeatureLocForDistance="TSS", output = "overlapping", maxgap = 2000, multiple = TRUE) %>% data.frame

xen.biomart <- biomaRt::getBM(attributes = 
                                c("ensembl_gene_id","external_gene_name","chromosome_name", 
                                  "start_position","end_position","description"
                                ),       
                              filters    = "",
                              values = "",
                              mart       = mart) 

t.brain.pca <- pca.brain %>%
  left_join(
    dplyr::select(xen.biomart, ensembl_gene_id, external_gene_name, description),
    by = c("feature" = "ensembl_gene_id"))

rlv.genes.brain <- c("grik2", "slc17a7", "grm4", "grm5", "grin1", "grm5", "grm1", "grm8", "grin2b", #glutamate signaling
                     "gabbr1", "gabbr2", "gabrb3", #GABA signaling
                     "dnmt3a", "hdac8", "mbd2", #methylation
                     "gh1", "irs2", "igfbp4", "igfbp5", #somatotropic
                     "trhr3", "trhde", "dio1", "\\btg\\b", #thyrotropic 
                     "esrrg",
                     "kiss2", "prkag2", "prkar1b", #GnRH signaling
                     "hsd17b12") #steroidogenesis


rlv.genes.table.brain <- t.brain.pca %>%
  filter(str_detect(external_gene_name, paste(rlv.genes.brain, collapse = "|")))

growth.genes.brain <- c(
                     "gh1", "irs2", "igfbp4", "igfbp5", #somatotropic
                     "trhr3", "trhde", "dio1", "\\btg\\b", #thyrotropic 
                     "esrrg")

growth.genes.table.brain <- t.brain.pca %>%
  filter(str_detect(external_gene_name, paste(growth.genes.brain, collapse = "|")))

t.growth.genes.table.brain <- growth.genes.table.brain %>% t %>% data.frame

  

pca.growth.genes.table.brain <- t.growth.genes.table.brain[grep(".Brain", rownames(t.growth.genes.table.brain)), ] %>% 
  rownames_to_column('ID') %>%
  mutate_at('ID', str_replace, "mcols.Brain_", "")

colnames(pca.growth.genes.table.brain) <- c("ID", t.growth.genes.table.brain["external_gene_name",])

growth.brain.table.pca <- merge(phenotype, pca.growth.genes.table.brain, by = "ID")

variables.pca <- growth.brain.table.pca[,c(3,4,8,13,17:31)] %>% mutate_if(is.character,as.numeric)

library(missMDA)
imp<-imputePCA(variables.pca)
imp.prop.variables<-imp$completeObs

library(factoextra)
library(FactoMineR)
pca.prop <- PCA(imp.prop.variables, scale.unit = TRUE, graph = FALSE)


##
biplot.genexpression<- fviz_pca_biplot(pca.prop, repel = TRUE,
                                       axes = c(1, 2),
                                       label = "var", habillage=as.factor(growth.brain.table.pca$treatment), addEllipses=FALSE, 
                                       col.var = "Gray31", # Variables color
                                       col.ind = c("Black"), 
                                       title = "PCA Growth Genes Brain"
)

biplot.genexpression       



######################################################################


PCA_phenotype <- function(meth.diff, perc.meth, rlv.genes, phenotype) {

mydiff <- getMethylDiff(meth.diff,
                                 difference=meth.cut,
                                 qvalue=qvalue.cut) %>% getData()

perc.meth.dm <- perc.meth[rownames(perc.meth) %in% rownames(mydiff),]

perc.meth.dm.merged <- cbind(mydiff, perc.meth.dm)

meth_pos.pca <- GRanges(
  seqnames = perc.meth.dm.merged$chr, 
  ranges = IRanges(start = perc.meth.dm.merged$start, end = perc.meth.dm.merged$end),
  mcols = perc.meth.dm.merged[, c(5:ncol(perc.meth.dm.merged))]
)

martfile <- c("./data/mart_xtropicalis_gene_ensembl.RData")

load(martfile)

xtro <- getAnnotation(mart, featureType = "TSS")

pca <- annotatePeakInBatch(meth_pos.pca, AnnotationData = xtro, FeatureLocForDistance="TSS", output = "overlapping", maxgap = 2000, multiple = TRUE) %>% data.frame

xen.biomart <- biomaRt::getBM(attributes = 
                                c("ensembl_gene_id","external_gene_name","chromosome_name", 
                                  "start_position","end_position","description"
                                ),       
                              filters    = "",
                              values = "",
                              mart       = mart) 

t.pca <- pca %>%
  left_join(
    dplyr::select(xen.biomart, ensembl_gene_id, external_gene_name, description),
    by = c("feature" = "ensembl_gene_id"))

rlv.genes.table <- t.pca %>%
  filter(str_detect(external_gene_name, paste(rlv.genes, collapse = "|")))

t.rlv.genes.table <- rlv.genes.table %>% t %>% data.frame



pca.rlv.genes.table <- t.rlv.genes.table[grep("", rownames(t.rlv.genes.table)), ] %>% 
  rownames_to_column('ID') %>%
  mutate_at('ID', str_replace, "mcols.", "")

colnames(pca.rlv.genes.table) <- c("ID", t.rlv.genes.table["external_gene_name",])

rlv.table.pca <- merge(phenotype, pca.rlv.genes.table, by = "ID")

variables.pca <- rlv.table.pca[,c(3:ncol(rlv.table.pca))] %>% mutate_if(is.character,as.numeric)

library(missMDA)
imp<-imputePCA(variables.pca)
imp.prop.variables<-imp$completeObs

library(factoextra)
library(FactoMineR)
pca.prop <- PCA(imp.prop.variables, scale.unit = TRUE, graph = FALSE)


##
biplot.genexpression<- fviz_pca_biplot(pca.prop, repel = TRUE,
                                       axes = c(1, 2),
                                       label = "var", habillage=as.factor(rlv.table.pca$treatment), addEllipses=FALSE, 
                                       col.var = "Gray31", # Variables color
                                       col.ind = c("Black"), 
                                       title = paste0("PCA", deparse(substitute(rlv.genes)))
)

biplot.genexpression     

}


phenotype <- read_excel("data/phenotype_data.xlsx")

brain.unite <- readRDS("./data/brain_united.rds")

pm.brain <- percMethylation(brain.unite) %>% 
  data.frame
colnames(pm.brain) <- str_replace(colnames(pm.brain), "Brain_", "")

brain.meth.dir <- c("./data/myDiff.brain.RData")
load(brain.meth.dir)

meth.cut <- 10
qvalue.cut <- 0.05


growth_genes <- c(
  "gh1", "irs2", "igfbp4", "igfbp5", #somatotropic
  "trhr3", "trhde", "dio1", "\\btg\\b", #thyrotropic 
  "esrrg")

size.phen <- phenotype %>%
  select(treatment, ID, body_weight,hindleg_length,glucose)



PCA_phenotype(myDiff.brain, pm.brain, growth_genes, size.phen)
