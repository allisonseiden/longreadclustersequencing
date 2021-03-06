# Felix Richter
# 11/18/2018
# Indel correlations with parental age
#################################################

setwd('D:/Dropbox/PhD/')
setwd('/Users/frichter/Dropbox (Personal)/PhD/')
setwd('/Users/felixrichter/Dropbox/PhD/')
options(stringsAsFactors=FALSE)

p = c('car', 'magrittr', 'extraDistr', 'ggplot2', 'tidyr', 'broom', 
      'readr', 'wesanderson', 'purrr', 'dplyr', 'stringi')
lapply(p, require, character.only = TRUE)

##
decode_indels = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/classified_indels.txt')
decode_dnvs = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/decode_DNMs.tsv')

decode_dnvs %>% head %>% as.data.frame
decode_dnvs$Phase_source %>% unique

## are the indels only for haplotype phasing? Yes
decode_dnvs %>% 
  mutate(indel = ifelse((stri_length(Ref) > 1) | (stri_length(Alt) > 1), 'indel', 'snv')) %>% 
  group_by(indel, Phase_source, Phase_combined) %>% tally

######################################################
# get IDs with 3-gen phasing
######################################################
three_gen_ids = decode_dnvs %>% 
  filter(Phase_source %in% c('three_generation', 'both_approaches')) %>% 
  select(Proband_nr) %>% unique %>% unlist %>% as.character
decode_other_ids = decode_dnvs %>% filter(!(Proband_nr %in% three_gen_ids)) %>% 
  select(Proband_nr) %>% unique %>% unlist %>% as.character
length(three_gen_ids) ## 225
length(decode_other_ids) ## 1323
decode_indels %>% filter(ID %in% three_gen_ids) %>% select(ID) %>% unique %>% dim
# 225 3-gen IDs but only 222 have indels
decode_indels %>% filter(ID %in% decode_other_ids) %>%
  select(ID) %>% unique %>% dim
# 1323 ILMN IDs but only 1308 have indels

# calculate fraction SNVs phased
fract_3gen_phased = decode_dnvs %>% 
  mutate(indel = ifelse((stri_length(Ref) > 1) | (stri_length(Alt) > 1), 'indel', 'snv')) %>% 
  filter(indel == 'snv') %>% 
  filter(Proband_nr %in% three_gen_ids) %>% 
  ## only keep 3-gen phased SNVs
  mutate(phase_3gen_only = Phase_source %in% c('three_generation', 'both_approaches')) %>% 
  # group_by(Phase_source, phase_3gen_only) %>% tally
  group_by(Proband_nr) %>% summarise(phased = sum(phase_3gen_only), total = n()) %>% ungroup %>% 
  mutate(fraction_snv_phased = phased/total) %>% 
  rename(ID = Proband_nr) %>% select(ID, fraction_snv_phased)
  
######################################################
# get an age dataframe
######################################################

age_df = decode_dnvs %>% 
    select(Proband_nr, Fathers_age_at_conception, Mothers_age_at_conception) %>% unique %>% 
    rename(father = Fathers_age_at_conception, mother = Mothers_age_at_conception, 
           ID = Proband_nr) %>% 
    gather(key = 'parent', value = 'parental_age', -ID) %>% unique %>% 
    mutate(parent = ifelse(parent == 'father', 'dad', 'mom')) %>% 
    mutate(cohort = ifelse(ID %in% three_gen_ids, 'decode_3gen', 'decode_ilmn'))

# write_tsv(age_df, 'longreadclustersequencing/literature/jonsson_decode_2017/age_df.txt')
age_df %>% group_by(cohort, parent) %>% tally
age_df_3gen = age_df %>% filter(cohort == 'decode_3gen')

# write_tsv(age_df_3gen, 'longreadclustersequencing/literature/jonsson_decode_2017/')
######################################
# Add in HR subtypes
######################################

decode_hr_grouping = decode_indels %>% 
  filter(Indel_Class == 'HR') %>% 
  mutate(hr_subtype = ifelse(Allele %in% c('A', 'T'), 'HR_AT', 'HR_GC')) %>% 
  select(ID:Alt, hr_subtype)

decode_indels %<>% left_join(decode_hr_grouping) 
# decode_indels %<>% mutate(Indel_Class = ifelse(is.na(hr_subtype), Indel_Class, hr_subtype))

######################################################
# join phasing info with indel class assignments
######################################################

indels_full_df = decode_dnvs %>% rename(Chrom = Chr, End = Pos_hg38, ID = Proband_nr) %>% 
  inner_join(decode_indels) %>% unique %>% 
  mutate(cohort = ifelse(ID %in% three_gen_ids, 'decode_3gen', 'decode_ilmn')) %>% 
  ## only keep the 3-gen IDs because none of the illumina indels seem to be phased.. 
  filter(ID %in% three_gen_ids)

indels_full_df$Phase_combined %>% unique
indels_full_df$Phase_source %>% unique
indels_full_df$ID %>% unique %>% length
indels_full_df %>% group_by(Indel_Class) %>% tally


indels_full_df
res_dir = 'longreadclustersequencing/phasing_analysis/results_phasing/'
write_tsv(indels_full_df, paste0(res_dir, 'indel_df_decode_2019_02_08.txt'))

########################################################
# new summary dataframe, accounting for IDs with
# 0 of a certain indel class
########################################################

indel_class_vec = c('CCC', 'non-CCC', 'HR') ## , 'HR_GC', 'HR_AT') # 

expand_df = expand.grid('ID' = age_df_3gen$ID, 'Indel_Class' = indel_class_vec,
                        # in_repeat = c('N', 'Y'),
                        'parent' = c('mom', 'dad')) %>% 
  full_join(age_df_3gen) %>% unique

## should be same
nrow(expand_df)
225*2*3
# 225*2*4*2

## summarise phasing results per proband
sum_indels = indels_full_df %>%
  rename(parent = Phase_combined) %>% 
  ## filter for non-repeats
  # filter(repFamily == '.') %>% 
  mutate(in_repeat = ifelse(repClass == '.', 'N', 'Y')) %>% 
  ## label the na 30/1241 as unphased
  # mutate(parent = ifelse(is.na(parent), 'unphased', parent)) %>% 
  # or just remove:
  filter(!is.na(parent)) %>% 
  ## summing
  mutate(parental_age = ifelse(parent == 'father',
                               Fathers_age_at_conception,
                               Mothers_age_at_conception)) %>%
  mutate(parent = ifelse(parent == 'father', 'dad', 'mom')) %>% 
  select(-Fathers_age_at_conception, -Mothers_age_at_conception) %>%
  group_by(cohort, ID, Indel_Class, parent, parental_age) %>% ## in_repeat
  summarise(Indel_ct = n()) %>% ungroup %>% 
  ## account for 0s
  full_join(expand_df) %>%
  mutate(Indel_ct = ifelse(is.na(Indel_ct), 0, Indel_ct)) 

# add an all indel category
decode_indel_all_df = sum_indels %>%
  group_by(ID, parent, parental_age, cohort) %>% ## in_repeat
  summarise(Indel_ct = sum(Indel_ct)) %>% ungroup %>% 
  mutate(Indel_Class = 'All')

sum_indels %<>% bind_rows(decode_indel_all_df)

## shoult add up to 225 per category
sum_indels %>% group_by(Indel_Class, parent) %>% tally

## join with 3-gen phased
decode_3gen_indels = sum_indels %>% inner_join(fract_3gen_phased)
write_tsv(decode_3gen_indels,
          'longreadclustersequencing/literature/jonsson_decode_2017/indel_cts_per_id_repeat_HRATGC_2019_02_06.txt')

####################
# Summary tables
####################

indels_full_df %>% 
  rename(parent = Phase_combined) %>% 
  mutate(parent = ifelse(is.na(parent), 'unphased', parent)) %>% 
  mutate(Indel_Class = ifelse(grepl('HR', Indel_Class), 'HR', Indel_Class)) %>% 
  mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% unique %>% 
  group_by(Indel_Class, ins_del, parent) %>% tally %>% 
  ungroup %>% 
  write_tsv('longreadclustersequencing/indel_analysis/summary_cts_decode.txt')


decode_dnvs %>% filter(!(Proband_nr %in% age_df_3gen$ID)) %>% 
  select(Proband_nr) %>% unique %>% nrow ## 1326

decode_indels %>% 
  filter(!(ID %in% age_df_3gen$ID)) %>% 
  mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% 
  mutate(Indel_Class = ifelse(grepl('HR', Indel_Class), 'HR', Indel_Class)) %>% unique %>% 
  group_by(Indel_Class, ins_del) %>% tally %>% ungroup %>% 
  write_tsv('longreadclustersequencing/indel_analysis/unphased_summary_cts_decode.txt')


#### plot the correlations
# p = sum_indels %>% 
#   ggplot(aes(x = parental_age, y = Indel_ct, color = parent)) +
#   geom_point() +
#   facet_wrap(~Indel_Class) +
#   scale_color_manual(values = c('#352846', '#E2AE38')) +
#   geom_smooth(method = 'lm') +
#   theme_classic()
# p
# 
# fit = sum_indels %>% 
#   # filter(Indel_ct < 10) %>% 
#   group_by(cohort, Indel_Class, parent) %>% 
#   do(model = lm(Indel_ct ~ parental_age, .))
# 
# fit %>% tidy(model) %>% filter(term != '(Intercept)')
# fit %>% glance(model)
# # fit %>% augment(model)
# 
# ## am I making a mistake?
# sum_indels$ID %>% unique %>% length
# sum_indels %>% unique %>% 
#   group_by(parent, Indel_Class) %>% tally
  
