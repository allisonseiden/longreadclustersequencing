# Felix Richter
# 1/22/2019
# Calculate the fraction phased
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
pb_indel_phased = paste0(res_dir, 'indels_df.txt') %>% read_tsv

ilmn_phased = paste0(res_dir, 'phasing_analysis_df_ilmn_2018_12_06.txt') %>% read_tsv

# NAs are from IDs that were not fully phased (remaining should be N=308)
ilmn_phased_ids = ilmn_phased %>% filter(!is.na(Mom)) %>% select(ID) %>% unique %>%
  unlist %>% as.character
pacbio_ids = c('1-00801', '1-01019', '1-03897', '1-04190', '1-04389',
               '1-04460', '1-04537', '1-05443', '1-05673', '1-05846')
pb_phased %>% filter(ID %in% pacbio_ids) %>% dim
pb_indel_phased %>% filter(ID %in% pacbio_ids) %>% dim

# Parental age
age_df = read_tsv('longreadclustersequencing/data/parental_age_at_conception_allPCGC_nonNA.txt')
# only have parental age for 305/208 Illumina trios
age_df_ilmn = age_df %>% filter(ID %in% ilmn_phased_ids) 
age_df_pb = age_df %>% filter(ID %in% pacbio_ids) 

## all have parental age data
age_df %>% filter(is.na(Maternal_age_at_conception) |
                    is.na(Paternal_age_at_conception)) %>% dim

##############################
# Calculate the fraction
# phased per sample
#############################

##############################
# for illumina
##############################

ilmn_fract_phased = ilmn_phased %>%
  filter(ID %in% age_df_ilmn$ID) %>% 
  filter(!is.na(Mom)) %>% 
  ## confirm there are none left in chrX and chrY
  # filter((Chrom %in% c('chrX', 'chrY')))
  mutate(indel = ifelse((stri_length(Ref) > 1) | (stri_length(Alt) > 1), 'indel', 'snv')) %>% 
  group_by(ID, indel) %>% 
  summarise(Mom = sum(Mom), Dad = sum(Dad), Unphased = sum(Unphased),
            total = Mom + Dad + Unphased) %>% ungroup %>% 
  mutate(fraction_phased = (total - Unphased)/total) %>% 
  select(ID, indel, fraction_phased) %>% 
  mutate(cohort = 'ilmn_305')

## account for the ID with 0 phased
full_ilmn_df = expand.grid('ID' = age_df_ilmn$ID, indel = c('snv', 'indel'))
ilmn_fract_phased %<>% full_join(full_ilmn_df) %>% 
  mutate(cohort = 'ilmn_305') %>% 
  mutate(fraction_phased = ifelse(is.na(fraction_phased), 0, fraction_phased))

##############################
# for PacBio
##############################

pb_fract_phased = pb_phased %>% 
  mutate(indel = ifelse((stri_length(Ref) > 1) | (stri_length(Alt) > 1), 'indel', 'snv')) %>% 
  mutate(IL_Unphased = ifelse((IL_Mom == 0) & (IL_Dad == 0) & (IL_Unphased == 0), 1, IL_Unphased)) %>% 
  ## add/remove ID depending on overall sum vs per ID 
  group_by(ID, indel) %>%
  summarise(PB_Mom = sum(PB_Mom), PB_Dad = sum(PB_Dad), PB_Unphased = sum(PB_Unphased),
            PB_total = PB_Mom + PB_Dad + PB_Unphased,
            IL_Mom = sum(IL_Mom), IL_Dad = sum(IL_Dad), IL_Unphased = sum(IL_Unphased),
            IL_total = IL_Mom + IL_Dad + IL_Unphased,
            true_total = n()) %>% ungroup %>% 
  mutate(fraction_phased_pb = (PB_total - PB_Unphased)/PB_total,
         fraction_phased_ilmn = (IL_total - IL_Unphased)/IL_total) %>% 
  select(ID, indel, fraction_phased_pb, fraction_phased_ilmn)

# for PacBio indels
pb_indel_sum = pb_indel_phased %>% 
  mutate(indel = 'indel') %>% 
  # filter((PB_Mom_Indel == 0) & (PB_Dad_Indel == 0) & (PB_Unphased_Indel == 0))
  group_by(ID, indel) %>%
  summarise(PB_Mom = sum(PB_Mom_Indel), PB_Dad = sum(PB_Dad_Indel),
            PB_Unphased = sum(PB_Unphased_Indel),
            PB_total = PB_Mom + PB_Dad + PB_Unphased,
            true_total = n()) %>% ungroup %>% 
  mutate(fraction_phased = (PB_total - PB_Unphased)/PB_total) %>% 
  select(ID, indel, fraction_phased) %>% 
  mutate(cohort = 'pb_10')

############################################################
# combine all data into a single dataframe
############################################################

ilmn_10_phased = pb_fract_phased %>% select(ID, indel, fraction_phased_ilmn) %>% 
  rename(fraction_phased = fraction_phased_ilmn) %>% 
  mutate(cohort = 'ilmn_10') 

fract_phased = pb_fract_phased %>% select(ID, indel, fraction_phased_pb) %>% 
  filter(indel != 'indel') %>% 
  rename(fraction_phased = fraction_phased_pb) %>% 
  mutate(cohort = 'pb_10') %>% 
  bind_rows(ilmn_10_phased, pb_indel_sum, ilmn_fract_phased)

fract_phased %>% group_by(cohort) %>% tally

write_tsv(fract_phased, 'longreadclustersequencing/phasing_analysis/fraction_phased/fraction_phased_2019_01_22.txt')

# weirdly, no correlation between fraction indel and snv phased
ilmn_fract_phased %>% spread(key = indel, value = fraction_phased) %>% 
  ggplot(aes(x = snv, y = indel)) + geom_point() + geom_smooth(method = 'lm')

## graph the distribution
p = fract_phased %>% 
  as.data.frame %>% 
  filter(cohort != 'ilmn_10') %>% 
  ggplot((aes(x = fraction_phased, fill = cohort))) +
  geom_histogram(alpha = 1, position = 'identity', bins = 40) +
  scale_fill_manual(values = c('grey60', 'red')) + ## 'red', blue
  facet_wrap(~indel, scales = 'free_y') +
  theme_classic()
p
ggsave('longreadclustersequencing/phasing_analysis/fraction_phased/fract_phased.png',
       p, width = 5, height = 2)

### some basic addition
phased_ilmn = (1217+3945+82+224)
unphased_ilmn = (16860 + 1229)
phased_ilmn/(phased_ilmn + unphased_ilmn)

### pb
phased_pb = (128+445+11+37)
unphased_pb = (105 + 10)
phased_pb/(phased_pb + unphased_pb)

