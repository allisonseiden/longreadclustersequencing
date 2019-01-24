
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

# local
p = c("magrittr", "stringi", "qvalue", "extraDistr", "purrr", "dplyr", "ggplot2", "tidyr", "readr") #
# minerva:
# library(extraDistr,lib="/hpc/users/richtf01/rLocalPackages/")
# p = c("magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "readr", "colorout") # 
lapply(p, require, character.only = TRUE)

dn = readRDS(file = "whole_genome/pcgcV2_trios/dn_allTSS_UpAndDownTSS_718mCHD_RefGeneUCSC_18_07_26.RDS")

parental_age_df = dn %>% filter(grepl('GMKF', wave)) %>% 
  select(Blinded.ID, Maternal.Age.at.Proband.Birth, Paternal.Age.at.Proband.Birth) %>% unique

parental_age_df %>% 
  filter(!is.na(Paternal.Age.at.Proband.Birth), !is.na(Maternal.Age.at.Proband.Birth)) %>% 
  mutate(Paternal_age_at_conception = Paternal.Age.at.Proband.Birth - 0.75,
         Maternal_age_at_conception = Maternal.Age.at.Proband.Birth - 0.75) %>% 
  rename(ID = Blinded.ID) %>% 
  select(ID, Maternal_age_at_conception, Paternal_age_at_conception) %>% 
  write_tsv("longreadclustersequencing/data/parental_age_at_conception_allPCGC_nonNA.txt")

