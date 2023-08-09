library(tidyverse)

qc <- read.csv("./data/multiqc_general_statistics.csv")

qc.bt <- qc %>%
  filter(str_detect(Sample.Name, "_T_|_B_")) %>%
  mutate(across(X..Aligned, parse_number))

mean(qc.bt$M.Seqs)
stats::sd(qc.bt$M.Seqs)

mean(qc.bt$X..Aligned)
stats::sd(qc.bt$X..Aligned)
