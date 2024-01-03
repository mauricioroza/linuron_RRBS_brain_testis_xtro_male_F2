
#' correlation_matrix
#' Creates a publication-ready / formatted correlation matrix, using `Hmisc::rcorr` in the backend.
#'
#' @param df dataframe; containing numeric and/or logical columns to calculate correlations for
#' @param type character; specifies the type of correlations to compute; gets passed to `Hmisc::rcorr`; options are `"pearson"` or `"spearman"`; defaults to `"pearson"`
#' @param digits integer/double; number of decimals to show in the correlation matrix; gets passed to `formatC`; defaults to `3`
#' @param decimal.mark character; which decimal.mark to use; gets passed to `formatC`; defaults to `.`
#' @param use character; which part of the correlation matrix to display; options are `"all"`, `"upper"`, `"lower"`; defaults to `"all"`
#' @param show_significance boolean; whether to add `*` to represent the significance levels for the correlations; defaults to `TRUE`
#' @param replace_diagonal boolean; whether to replace the correlations on the diagonal; defaults to `FALSE`
#' @param replacement character; what to replace the diagonal and/or upper/lower triangles with; defaults to `""` (empty string)
#'
#' @return a correlation matrix
#' @export
#'
#' @examples
#' `correlation_matrix(iris)`
#' `correlation_matrix(mtcars)`
correlation_matrix <- function(df, 
                               type = "pearson",
                               digits = 3, 
                               decimal.mark = ".",
                               use = "all", 
                               show_significance = TRUE, 
                               replace_diagonal = FALSE, 
                               replacement = ""){
  
  # check arguments
  stopifnot({
    is.numeric(digits)
    digits >= 0
    use %in% c("all", "upper", "lower")
    is.logical(replace_diagonal)
    is.logical(show_significance)
    is.character(replacement)
  })
  # we need the Hmisc package for this
  require(Hmisc)
  
  # retain only numeric and boolean columns
  isNumericOrBoolean = vapply(df, function(x) is.numeric(x) | is.logical(x), logical(1))
  if (sum(!isNumericOrBoolean) > 0) {
    cat('Dropping non-numeric/-boolean column(s):', paste(names(isNumericOrBoolean)[!isNumericOrBoolean], collapse = ', '), '\n\n')
  }
  df = df[isNumericOrBoolean]
  
  # transform input data frame to matrix
  x <- as.matrix(df)
  
  # run correlation analysis using Hmisc package
  correlation_matrix <- Hmisc::rcorr(x, type = )
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  # transform correlations to specific character format
  Rformatted = formatC(R, format = 'f', digits = digits, decimal.mark = decimal.mark)

  # add significance levels if desired
  if (show_significance) {
    # define notions for significance levels; spacing is important.
    stars <- ifelse(is.na(p), "   ", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   "))))
    Rformatted = paste0(Rformatted, stars)
  }
  # build a new matrix that includes the formatted correlations and their significance stars
  Rnew <- matrix(Rformatted, ncol = ncol(x))
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep =" ")
  
  # replace undesired values
  if (use == 'upper') {
    Rnew[lower.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (use == 'lower') {
    Rnew[upper.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (replace_diagonal) {
    diag(Rnew) <- replacement
  }
  
  return(Rnew)
}



corr_phenotype <- function(meth.diff, perc.meth, rlv.genes, phenotype) {
  
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
    filter(str_detect(external_gene_name, paste(rlv.genes, collapse = "|"))) %>%
    arrange(start)%>%
    arrange(factor(external_gene_name, levels = rlv.genes))
  
  t.rlv.genes.table <- rlv.genes.table %>% t %>% data.frame
  
  
  
  pca.rlv.genes.table <- t.rlv.genes.table[grep("", rownames(t.rlv.genes.table)), ] %>% 
    rownames_to_column('ID') %>%
    mutate_at('ID', str_replace, "mcols.", "")
  
  colnames(pca.rlv.genes.table) <- c("ID", t.rlv.genes.table["external_gene_name",])
  
  rlv.table.pca <- merge(phenotype, pca.rlv.genes.table, by = "ID")
  
  variables.pca <- rlv.table.pca[,c(3:ncol(rlv.table.pca))] %>% mutate_if(is.character,as.numeric)
  
  variables.pca
}

brain_genes <- c("grik2", "slc17a7", "grm4", "grm5", "grin1", "grm1", "grm8", "grin2b", #glutamate signaling
                                    "gabbr1", "gabbr2", "gabrb3", #GABA signaling
                                    "gh1", "irs2", "igfbp4", "igfbp5", #somatotropic
                                    "trhr3", "trhde", "dio1", "\\btg\\b", #thyrotropic 
                                    "kiss2", "prkag2", "prkar1b", #GnRH signaling
                                    "hsd17b12", "esrrg", #steroidogenesis
                                    "dnmt3a", "hdac8", "mbd2" #methylation
                 )

phen <- phenotype %>%
  dplyr::select(treatment, ID, body_weight, svl_length, hindleg_length, germ_cells_nests, fertility_rate)

df_brain <- corr_phenotype(myDiff.brain, pm.brain, brain_genes, phen)


corr_brain <- correlation_matrix(df_brain, type = "spearman", digits = 2)

corr_brain <- corr_brain[-c(1:(ncol(phen)-2)), -c((ncol(phen)-1):ncol(corr_brain))] %>% data.frame %>%
  rownames_to_column(var = "Gene")


#######
# Testis

testis_genes <- c("piwil1", "mael", "spo11", "\\bddx4\\b", #spermatogenesis, gonadal development
                  "dnmt3a", "ep300", "elp3", "kat5", "kat14", #histone acetyltransferase, methylation
                  "hsd17b12", "esrrg", "lhcgr")

df_testis <- corr_phenotype(myDiff.testis, pm.testis, testis_genes, phen)

corr_testis <- correlation_matrix(df_testis, type = "spearman", digits = 2)

corr_testis <- corr_testis[-c(1:(ncol(phen)-2)), -c((ncol(phen)-1):ncol(corr_testis))] %>% data.frame %>%
  rownames_to_column(var = "Gene")

corr_dfs <- list("Brain" = corr_brain, "Testis" = corr_testis)

library(writexl)
write_xlsx(corr_dfs, path = "./supplementary_material_tables/correlations.xlsx")

