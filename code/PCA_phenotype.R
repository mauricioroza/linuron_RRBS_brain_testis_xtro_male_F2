library(tidyverse)
library(readxl)
library(ChIPpeakAnno)
library(methylKit)

######################################################################


PCA_phenotype <- function(meth.diff, perc.meth, rlv.genes, phenotype, plot_title) {

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
                                         title = plot_title
  )
  
  biplot.genexpression     

}


phenotype <- read_excel("data/phenotype_data.xlsx") %>%
  mutate(fertility_rate = fertilized_eggs/total_eggs)

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
  "kiss2", "prkag2", "prkar1b", #GnRH signaling
  "hsd17b12" #steroidogenesis
  )

size.phen <- phenotype %>%
  dplyr::select(treatment, ID, body_weight,hindleg_length,glucose)



PCA_phenotype(myDiff.brain, pm.brain, growth_genes, size.phen, "PCA Growth Genes Brain")

reproduction_genes <- c("kiss2", "prkag2", "prkar1b", #GnRH signaling
                        "hsd17b12") #steroidogenesis

rep_phen <- phenotype %>%
  dplyr::select(treatment, ID,glucose, germ_cells_nests)

PCA_phenotype(myDiff.brain, pm.brain, reproduction_genes, rep_phen, "PCA Reproduction Genes Brain")

###########################################################
#Testis genes

testis.unite <- readRDS("./data/testis_unite.rds")

pm.testis <- percMethylation(testis.unite) %>% 
  data.frame
colnames(pm.testis) <- str_replace(colnames(pm.testis), "testis_", "")

testis.meth.dir <- c("./data/myDiff.testis.RData")
load(testis.meth.dir)

testis_genes <- c("lhcgr", "hsd17b12", "esrrg",
                  "piwil1", "mael", "spo11", "\\bddx4\\b" #spermatogenesis, gonadal development
)

testis_phen <- phenotype %>%
  dplyr::select(treatment, ID, germ_cells_nests, fertility_rate)

PCA_phenotype(myDiff.testis, pm.testis, testis_genes, testis_phen, "PCA Reproduction Genes Testis")



