library(rstatix)
library(writexl)

phenotype <- read_excel("data/phenotype_data.xlsx") %>%
  mutate(fertility_rate = fertilized_eggs/total_eggs,
         BMI = (body_weight)/(svl_length^2)*1000,
         hll_svl = hindleg_length/svl_length)

phenotype <- phenotype %>%
  dplyr::select(treatment, ID, body_weight, svl_length, BMI, hindleg_length, hll_svl, glucose, germ_cells_nests)


tests <- lapply(phenotype[,c(3:ncol(phenotype))], function (x) {wilcox.test(x ~ treatment, data = phenotype, exact = FALSE)})

# tests <- lapply(phenotype[,c(3:ncol(phenotype))], function (x) {lm(x ~ treatment, data = phenotype)})
# 
# summ <- lapply(tests, function (x) {summary(x)})
# summ

summary <- phenotype %>% 
  group_by(treatment)%>%
  get_summary_stats(type = "mean_se")

summ_ctrl <- subset(summary, summary$treatment == "Control")
summ_lin <- subset(summary, summary$treatment == "Linuron_H")


stars_sig_table <- function(tests_list, summ_list,
                               
                               digits = 2, 
                               decimal.mark = "."){
  
  # check arguments
  stopifnot({
    is.numeric(digits)
    digits >= 0
  })

  p_values <- sapply(tests_list, function(p) p$p.value)
  p <- data.frame(p.value = p_values, row.names = names(tests_list))
  
  summ_list_mean <- summ_list$mean
  
  # transform correlations to specific character format
  Rformatted = formatC(summ_list_mean, format = 'f', digits = digits, decimal.mark = decimal.mark)
  
  # add significance levels if desired
  if (show_significance) {
    # define notions for significance levels; spacing is important.
    stars <- ifelse(is.na(p), "   ", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   "))))
    Rformatted = paste0(Rformatted, stars)
  }
  
  new_table <- summ_list %>%
    mutate(mean = Rformatted)

  return(new_table)
}

summ_lin_stars <- stars_sig_table(tests, summ_lin)

sig_tables <- list(Control = summ_ctrl, Linuron = summ_lin_stars)

write_xlsx(sig_tables, path = "./supplementary_material_tables/phenotype_table.xlsx")

