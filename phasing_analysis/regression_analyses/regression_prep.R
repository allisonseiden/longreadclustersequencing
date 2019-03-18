# Felix Richter
# 1/22/2019
# Prepare dataframe 
# for indel correlations with parental age
#################################################

setwd('D:/Dropbox/PhD/')
setwd('/Users/frichter/Dropbox (Personal)/PhD/')
setwd('/Users/felixrichter/Dropbox/PhD/')
options(stringsAsFactors=FALSE)

p = c('magrittr', 'extraDistr', 'broom', 'car',
      'purrr', 'dplyr', 'ggplot2', 'tidyr', 'readr', 'wesanderson', 'stringi')
lapply(p, require, character.only = TRUE)

###########################
# Original phasing results
###########################

res_dir = 'longreadclustersequencing/phasing_analysis/results_phasing/'
pb_phased = paste0(res_dir, 'phasing_analysis_df.txt') %>% read_tsv
ilmn_phased = paste0(res_dir, 'phasing_analysis_df_ilmn_2018_12_06.txt') %>% read_tsv

# get the list of phased IDs
# NAs are from IDs that were not fully phased (remaining should be N=308)
ilmn_phased_ids = ilmn_phased %>% filter(!is.na(Mom)) %>% select(ID) %>% unique %>%
  unlist %>% as.character
pacbio_ids = c('1-00801', '1-01019', '1-03897', '1-04190', '1-04389',
               '1-04460', '1-04537', '1-05443', '1-05673', '1-05846')
pb_phased %>% filter(ID %in% pacbio_ids) %>% dim
# pb_indel_phased %>% filter(ID %in% pacbio_ids) %>% dim

# cbind('phased_ids' = c(ilmn_phased_ids, pacbio_ids)) %>% as.data.frame %>% 
#   arrange(phased_ids) %>% 
#   write_tsv('longreadclustersequencing/data/phased_ids_pcgc.txt')

###############
# Parental age
###############

age_df = read_tsv('longreadclustersequencing/data/parental_age_at_conception_allPCGC_nonNA.txt')
# only have parental age for 305/208 Illumina trios
age_df_ilmn = age_df %>% filter(ID %in% ilmn_phased_ids) 
age_df_pb = age_df %>% filter(ID %in% pacbio_ids) 

## all have parental age data
age_df %>% filter(is.na(Maternal_age_at_conception) |
                    is.na(Paternal_age_at_conception)) %>% dim

##############################
# Load ancestry PCs
#############################

## only have ancestry for 23 samples (those with RNAseq)
# pca_evec = read_tsv("whole_genome/ancestry/pca_results_wgs_s116_t350_t416/wgs_s116_t350_t416_1kg_merged.Unsupervised.pca.evec.txt",
#                     col_names = F)
# names(pca_evec) = c("rm", "Blinded.ID", paste0("PC", seq_along(1:21)))
# pca_evec %>% head
# pca_values = pca_evec %>%
#   select(-rm, -PC21) %>% 
#   filter(grepl("1-", Blinded.ID)) %>% 
#   separate(Blinded.ID, "ID", sep = ":", extra = "drop") %>% 
#   filter(ID %in% age_df_ilmn$ID)

######################################
# Load and clean joined indel results
# and join with parental age
######################################

pb_indel_class_df = read_tsv('longreadclustersequencing/indel_analysis/all_indel_info.txt') %>% 
  select(-X1)
pb_indel_class_df$ID %>% unique %>% length

ilmn_indel_df = paste0(res_dir, 'indels_df_ilmn_2018_12_08.txt') %>% read_tsv
# note that 1 Illumina ID had 0 indels

# Join with age
# why are there NA indels? Most are IDs that did not complete phasing
ilmn_indel_df %<>% right_join(age_df_ilmn)
pb_indel_class_df %<>% inner_join(age_df_pb)

## the NA chrom is because of the ID without any indels (so assign Mom, Dad, Unphased to 0)
## the remaining phased NA are chrX and chrY so remove all chrX and chrY
ilmn_indel_df %>% filter(is.na(Mom)) %>% group_by(Chrom) %>% tally

ilmn_indel_df %<>%
  ## remove chrX and chrY (even if some are phased, just in case)
  filter(!(Chrom %in% c('chrX', 'chrY'))) %>%
  ## remove the two NA indels (chr1 and chr19), keeping the ID without any indels
  filter((!is.na(Mom)) | is.na(Chrom)) %>% 
  ## set the ID without indels as 0 for key variables
  mutate_at(vars(Mom:Unphased), function(x) ifelse(is.na(x), 0, x))

######################################
# Create pacbio indel class column
######################################

pb_indel_class_df %<>% 
  mutate(Indel_Class = ifelse(HR == 1, 'HR', NA)) %>% 
  mutate(Indel_Class = ifelse(CCC == 1, 'CCC', Indel_Class)) %>% 
  mutate(Indel_Class = ifelse(`non-CCC` == 1, 'non-CCC', Indel_Class))
  
######################################
# Add in HR subtypes
######################################

ilmn_hr_grouping = ilmn_indel_df %>% 
  filter(Indel_Class == 'HR') %>% 
  mutate(hr_subtype = ifelse(Allele %in% c('A', 'T'), 'HR_AT', 'HR_GC')) %>% 
  select(ID:Alt, hr_subtype)

# ilmn_indel_df %<>% left_join(ilmn_hr_grouping) 

pb_hr_grouping = pb_indel_class_df %>% 
  filter(HR == 1) %>% 
  mutate(hr_subtype = ifelse(Allele %in% c('A', 'T'), 'HR_AT', 'HR_GC')) %>% 
  select(ID:Location, Ref, Alt, hr_subtype)

# pb_indel_class_df %<>% left_join(pb_hr_grouping) 

# ilmn_indel_df %<>% mutate(Indel_Class = ifelse(is.na(hr_subtype), Indel_Class, hr_subtype))
# pb_indel_class_df %<>% mutate(Indel_Class = ifelse(is.na(hr_subtype), Indel_Class, hr_subtype))

write_tsv(ilmn_indel_df, paste0(res_dir, 'indel_df_ilmn_2019_02_07.txt'))
write_tsv(pb_indel_class_df, paste0(res_dir, 'indel_df_pb_2019_02_07.txt'))


######################################
# Sum the number of variants in 
# each variant class per ID
######################################

ilmn_indel_sum = ilmn_indel_df %>% 
  mutate(in_repeat = ifelse(repClass == '.', 'N', 'Y')) %>% 
  # mutate(repClassSuper = ifelse(repClass %in% c('LINE', 'SINE', 'LTR'), repClass, 'Other')) %>% 
  # mutate(repClassSuper = ifelse(in_repeat == 'Y', repClassSuper, '.')) %>% 
  ## only keep phased indels
  filter((Mom == 1) | (Dad == 1)) %>% 
  mutate(parent = ifelse(Mom == 1, 'mom', 'dad')) %>% 
  mutate(parental_age = ifelse(parent == 'dad',
                               Paternal_age_at_conception,
                               Maternal_age_at_conception)) %>%
  group_by(ID, Indel_Class, parent, parental_age, in_repeat) %>% ## 
  summarise(Indel_ct = n()) %>%
  ungroup

pb_indel_sum = pb_indel_class_df %>% 
  ## account for repeats
  mutate(in_repeat = ifelse(repClass == '.', 'N', 'Y')) %>% 
  # mutate(repClassSuper = ifelse(repClass %in% c('LINE', 'SINE', 'LTR'), repClass, 'Other')) %>% 
  # mutate(repClassSuper = ifelse(in_repeat == 'Y', repClassSuper, '.')) %>% 
  ## assign parents
  filter((PB_Mom_Indel == 1) | (PB_Dad_Indel == 1)) %>% 
  mutate(parent = ifelse(PB_Mom_Indel == 1, 'mom', 'dad')) %>% 
  mutate(parental_age = ifelse(parent == 'dad',
                               Paternal_age_at_conception,
                               Maternal_age_at_conception)) %>%
  group_by(ID, Indel_Class, parent, parental_age, in_repeat) %>% 
  summarise(Indel_ct = n()) %>%
  ungroup

indel_sum = bind_rows('ilmn_305' = ilmn_indel_sum, 'pb_10' = pb_indel_sum, .id = 'cohort')

######################################
# indel expand DF (get full dataframe accounting for 0 observed)
######################################

indel_class_vec = c('CCC', 'non-CCC', 'HR') ## 'HR_GC', 'HR_AT') #
# need dataframe with columns:
# ID, parent, Indel_Class, Indel_ct, parental_age, fraction_snvs_phased, ancestry
# nrows = length(ID) x 2 parents x 3 Indel_Class x 2 for repeat
305*2*4*2
id_class_df = expand.grid('ID' = age_df_ilmn$ID, 'Indel_Class' = indel_class_vec,
                          in_repeat = c('N', 'Y'),
                          'parent' = c('mom', 'dad')) %>% full_join(age_df_ilmn) %>% 
  mutate(parental_age = ifelse(parent == 'dad',
                               Paternal_age_at_conception,
                               Maternal_age_at_conception)) %>% 
  select(-Paternal_age_at_conception, -Maternal_age_at_conception)
nrow(id_class_df)

pb_id_class_df = expand.grid('ID' = age_df_pb$ID, 'Indel_Class' = indel_class_vec,
                             in_repeat = c('N', 'Y'),
                             'parent' = c('mom', 'dad')) %>% full_join(age_df_pb) %>% 
  mutate(parental_age = ifelse(parent == 'dad',
                               Paternal_age_at_conception,
                               Maternal_age_at_conception)) %>% 
  select(-Paternal_age_at_conception, -Maternal_age_at_conception)

nrow(pb_id_class_df)

id_class_df = bind_rows('ilmn_305' = id_class_df, 'pb_10' = pb_id_class_df, .id = 'cohort')
nrow(id_class_df)

# get the full dataframe
indel_full_ct = indel_sum %>% full_join(id_class_df) %>%
  mutate(Indel_ct = ifelse(is.na(Indel_ct), 0, Indel_ct)) %>% 
  arrange(parent, Indel_Class, ID)

# add an 'All' category
indel_all_df = indel_full_ct %>% group_by(ID, parent, parental_age, cohort, in_repeat) %>%
  summarise(Indel_ct = sum(Indel_ct)) %>% ungroup %>% 
  mutate(Indel_Class = 'All')

indel_full_ct %<>% bind_rows(indel_all_df)

indel_full_ct %>% group_by(cohort) %>% tally

write_tsv(indel_full_ct, 'longreadclustersequencing/indel_analysis/counts_per_id_ilmn_pb_repeat_HR_AT_GC_2019_02_06.txt')

ilmn_indel_df %>%
  mutate(Indel_Class = ifelse(grepl('HR', Indel_Class), 'HR', Indel_Class)) %>% 
  mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% 
  group_by(Indel_Class, ins_del, Mom, Dad, Unphased) %>% tally %>% 
  ungroup %>% 
  write_tsv('longreadclustersequencing/indel_analysis/summary_cts_ilmn_305.txt')
  
pb_indel_class_df %>% 
  mutate(Indel_Class = ifelse(grepl('HR', Indel_Class), 'HR', Indel_Class)) %>% 
  mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% 
  group_by(Indel_Class, ins_del, PB_Mom_Indel, PB_Dad_Indel, PB_Unphased_Indel) %>% tally %>% 
  ungroup %>% 
  write_tsv('longreadclustersequencing/indel_analysis/summary_cts_pb_10.txt')

############################################
# compare to the python dataframes
############################################
# indel_df_dad = read_tsv('longreadclustersequencing/phasing_analysis/indel_df_troubleshooting/indels_dad.txt')
# indel_df_mom = read_tsv('longreadclustersequencing/phasing_analysis/indel_df_troubleshooting/indels_mom.txt')
# 
# indel_df_dad %>% select(Paternal_age_at_conception:All) %>% as.matrix %>% cor
# cor.test(indel_df_dad$Paternal_age_at_conception, indel_df_dad$CCC, use = 'complete.obs')
# 
# indel_df_dad_subset = indel_df_dad %>% 
#   select(ID, Paternal_age_at_conception, CCC) %>% 
#   rename(dad_age_python = Paternal_age_at_conception, CCC_python = CCC) %>% unique
# 
# py_r_ccc_dad_joined = ilmn_indel_full_ct %>%
#   filter(parent == 'dad', Indel_Class == 'CCC') %>% 
#   full_join(indel_df_dad_subset)
# 
# py_r_ccc_dad_joined %>% filter(is.na(Indel_ct)) %>% as.data.frame
# 
# py_r_ccc_dad_joined %>% 
#   filter(!is.na(Indel_ct)) %>% 
#   summarise(r = cor(Indel_ct, parental_age, method = 'pearson', use = 'complete.obs'),
#             r_py_ccc = cor(CCC_python, parental_age, method = 'pearson', use = 'complete.obs'),
#             r_py_age = cor(Indel_ct, dad_age_python, method = 'pearson', use = 'complete.obs'),
#             r_py = cor(CCC_python, dad_age_python, method = 'pearson', use = 'complete.obs'),
#             p = cor.test(Indel_ct, parental_age)$p.value,
#             p_py = cor.test(CCC_python, dad_age_python, method = 'pearson', use = 'complete.obs')$p.value,
#             Indel_ct = sum(Indel_ct, na.rm = T))
# 
# py_r_ccc_dad_joined %>% 
#   filter(Indel_ct == CCC_python)
#   # filter(parental_age < (dad_age_python - 0.1))
#   filter(parental_age > (dad_age_python + 0.1))
# 


# ######################################
# # Investigate the proportion of DNVs
# # in specific repeats
# ######################################
# ilmn_indels_per_repClass = ilmn_indel_df %>% 
#   filter(repClass != '.') %>%
#   filter((Mom == 1) | (Dad == 1)) %>% 
#   mutate(parent = ifelse(Mom == 1, 'mom', 'dad')) %>% 
#   mutate(parental_age = ifelse(parent == 'dad',
#                                Paternal_age_at_conception,
#                                Maternal_age_at_conception)) %>%
#   mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% 
#   # filter(Indel_Class == 'CCC') %>% #head %>% as.data.frame
#   mutate(repClassSuper = ifelse(repClass %in% c('LINE', 'SINE', 'LTR'), repClass, 'Other')) %>% 
#   # filter(parent == 'mom') %>% 
#   # filter(repClass == 'LTR') %>% as.data.frame
#   ### ins_del, 
#   group_by(Indel_Class, parent, ins_del, repClassSuper) %>%
#   summarise(Indel_Ct = n()) %>% ungroup
