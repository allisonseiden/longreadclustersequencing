# Felix Richter
# 2/07/2019
# Replication timing
#################################################

setwd('D:/Dropbox/PhD/')
setwd('/Users/frichter/Dropbox (Personal)/PhD/')
setwd('/Users/felixrichter/Dropbox/PhD/')
options(stringsAsFactors=FALSE)

p = c('magrittr', 'extraDistr', 'broom', 'car',
      'purrr', 'dplyr', 'ggplot2', 'tidyr', 'readr', 'wesanderson', 'stringi')
lapply(p, require, character.only = TRUE)

####################################################
# load and clean replication timing data
####################################################

rep_timing_f_list = list.files('longreadclustersequencing/data/repliseq_anno', 'ENC.*.bed',
                               full.names = T)
rep_timing_f_list = list.files('longreadclustersequencing/literature/jonsson_decode_2017/repliseq_anno',
                               'ENC.*.bed', full.names = T)
names(rep_timing_f_list) = rep_timing_f_list %>% gsub('.*sorted_hg19_|.bed','', .)
col_name_list = c('Chrom_hg19', 'Start_hg19', 'End_hg19', 'var_id', 'var_len', 'RT')
rt_df = map_df(rep_timing_f_list, read_tsv, col_names = col_name_list, .id = 'RT_data')

phase_map = read_tsv('longreadclustersequencing/data/repliseq_anno/sample_phase_map.txt',
                     col_names = c('RT_data', 'Phase'))

max_rt_per_var = rt_df %>%
  inner_join(phase_map) %>% 
  group_by(RT_data) %>% 
  mutate(RT_z = (RT - mean(RT))/sd(RT)) %>% ungroup %>% 
  ## calculate difference to the next closest
  arrange(var_id, desc(RT_z)) %>% 
  mutate(diff_RT = c(diff(RT_z), NA)) %>% 
  ## Pick the maximum RT per variant ID
  group_by(var_id) %>% 
  filter(RT_z == max(RT_z)) %>% 
  ungroup

## only pick values that are clearly higher in one vs the next best (e.g., 1 standard deviation)
max_rt_per_var %<>% filter(abs(diff_RT) >= pnorm(2))

  ## check out the min/max of the RT z-scores
max_rt_per_var %>% ggplot(aes(x = RT_z)) + geom_histogram()
max_rt_per_var %>% filter(RT_z < 0) %>% 
  select(var_id, RT_z) %>% inner_join(rt_df) %>% inner_join(phase_map)

## exclude the two variants with 0 reads in all Repliseq experiments
max_rt_per_var %<>% filter(RT_z >= 0) 

max_rt_per_var %>% 
  # group_by(var_id) %>% tally %>% ungroup %>% filter(n > 1)
  group_by(Phase) %>% tally
  
# rt_wide = rt_df %>% 
#   select(-var_len) %>% 
#   spread(key = RT_data, value = RT)
# 
# ## confirm no NAs
# is.na(rt_wide) %>% sum
# 
# ## scale each RT column
# rt_wide_z = rt_df %>% 
#   select(-var_len) %>% 
#   group_by(RT_data) %>% 
#   mutate(RT_z = (RT - mean(RT))/sd(RT)) %>% ungroup %>% select(-RT) %>% 
#   spread(key = RT_data, value = RT_z)
# 
# rt_wide %>% head %>% as.data.frame
# rt_wide_z %>% head %>% as.data.frame
# 
# write_tsv(rt_wide, 'longreadclustersequencing/data/repliseq_anno/rt_wide_decode.txt')
write_tsv(max_rt_per_var, 'longreadclustersequencing/data/repliseq_anno/rt_assigned_Diff2sd_decode.txt')
# write_tsv(max_rt_per_var, 'longreadclustersequencing/data/repliseq_anno/rt_assigned_Diff2sd.txt')

## High in a certain RT == recently replicated at that point in time
## G1 → S → G2

max_rt_per_var = read_tsv('longreadclustersequencing/data/repliseq_anno/rt_assigned_hiDiff.txt')
max_rt_per_var_decode = read_tsv('longreadclustersequencing/data/repliseq_anno/rt_assigned_hiDiff_decode.txt')


####################################################
# Compare replication data between classes
# and between young vs old fathers w/in
# same quintile
####################################################

res_dir = 'longreadclustersequencing/phasing_analysis/results_phasing/'
ilmn_indel_df = read_tsv(paste0(res_dir, 'indel_df_ilmn_2019_02_07.txt'))
pb_indel_class_df = read_tsv(paste0(res_dir, 'indel_df_pb_2019_02_07.txt'))
decode_indel_df = read_tsv(paste0(res_dir, 'indel_df_decode_2019_02_08.txt'))

ilmn_rt_df = ilmn_indel_df %>% 
  unite(col = 'var_id', Chrom, Start, End, Ref, Alt, ID, sep = '.', remove = F) %>% 
  inner_join(max_rt_per_var) %>% unique %>% 
  mutate(parent = ifelse(Mom == 1, 'Mother', 'Father')) %>%
  # filter(Unphased == 0) %>%
  mutate(parent = ifelse(Unphased == 1, 'Unphased', parent)) %>% 
  mutate(parent_age = ifelse(parent == 'Mother', Maternal_age_at_conception, Paternal_age_at_conception)) %>% 
  mutate(parent_age = ifelse(parent == 'Unphased', NA, parent_age)) %>% 
  mutate(cohort = 'ilmn_305') %>% 
  select(var_id, ID, parent_age, Indel_Class, cohort, parent, Phase, RT_z)

pb_rt_df = pb_indel_class_df %>% 
  mutate(Start = Location - 1, End = Location) %>% 
  unite(col = 'var_id', Chrom, Start, End, Ref, Alt, ID, sep = '.', remove = F) %>% 
  inner_join(max_rt_per_var) %>% unique %>% 
  mutate(parent = ifelse(PB_Mom_Indel == 1, 'Mother', 'Father')) %>%
  # filter(Unphased == 0) %>%
  mutate(parent = ifelse(PB_Unphased_Indel == 1, 'Unphased', parent)) %>% 
  mutate(parent_age = ifelse(parent == 'Mother', Maternal_age_at_conception, Paternal_age_at_conception)) %>% 
  mutate(parent_age = ifelse(parent == 'Unphased', NA, parent_age)) %>% 
  mutate(cohort = 'pb_10') %>% 
  select(var_id, ID, parent_age, Indel_Class, cohort, parent, Phase, RT_z)

decode_rt_df  =decode_indel_df %>%
  unite(col = 'var_id', Chrom, Start, End, Ref, Alt, ID, sep = '.', remove = F) %>% 
  inner_join(max_rt_per_var_decode) %>% unique %>% 
  mutate(parent = ifelse(Phase_combined == 'mother', 'Mother', 'Father')) %>%
  mutate(parent = ifelse(is.na(Phase_combined), 'Unphased', parent)) %>%
  # group_by(Phase_combined, parent) %>% tally
  mutate(parent_age = ifelse(parent == 'Mother', Mothers_age_at_conception, Fathers_age_at_conception)) %>% 
  mutate(parent_age = ifelse(parent == 'Unphased', NA, parent_age)) %>% 
  select(var_id, ID, parent_age, Indel_Class, cohort, parent, Phase, RT_z)

#######################################################
# Question: are CCCs earlier in RT compared to others?
#######################################################

GetIndelsPerPhase = function(rt_df) {
  indels_per_phase = rt_df %>% 
    group_by(Indel_Class, parent) %>% #
    mutate(n_class = n()) %>% ungroup %>% 
    group_by(Indel_Class, Phase, parent) %>% ## , parent
    summarise(frac_indels_in_phase = n()/unique(n_class)) %>% ungroup
  
  indels_all_per_phase = rt_df %>% 
    group_by(parent) %>% #
    mutate(n_class = n()) %>% ungroup %>% 
    group_by(Phase, parent) %>% ## , parent
    summarise(frac_indels_in_phase = n()/unique(n_class)) %>% ungroup %>% mutate(Indel_Class = 'All')
  indels_per_phase %<>% bind_rows(., indels_all_per_phase)
  
  return(indels_per_phase)
}

ilmn_indels_per_phase  = GetIndelsPerPhase(ilmn_rt_df) 
pb_indels_per_phase  = GetIndelsPerPhase(pb_rt_df) 
decode_indels_per_phase  = GetIndelsPerPhase(decode_rt_df) 

indels_per_phase = bind_rows('ilmn' = ilmn_indels_per_phase,
                             'pb' = pb_indels_per_phase,
                             'decode' = decode_indels_per_phase,
                             .id = 'cohort')

phase_order = unique(indels_per_phase$Phase)[c(1, 3:6, 2)]
phase_order
class_order = c('HR', 'CCC', 'non-CCC', 'All')

## account for 0s
expand_df = expand.grid('Phase' = phase_order,
                        'parent' = c('Father', 'Mother', 'Unphased'),
                        cohort = unique(indels_per_phase$cohort),
                        'Indel_Class' = class_order)
indels_per_phase %<>% full_join(expand_df) %>% 
  mutate(frac_indels_in_phase = ifelse(is.na(frac_indels_in_phase), 0, frac_indels_in_phase))

## Question: are CCCs earlier in RT compared to others?
p = indels_per_phase %>%
  mutate(Indel_Class = factor(Indel_Class, levels = class_order)) %>% 
  mutate(Phase = factor(Phase, levels = phase_order)) %>% 
  ggplot(aes(x = Phase, y = frac_indels_in_phase, color = Indel_Class, fill = Indel_Class,
             group = Indel_Class)) +
  # geom_col(position = 'dodge2') +
  # geom_line() +
  geom_point(stat='summary', fun.y=sum) +
  stat_summary(fun.y=sum, geom="line") + 
  ## adding a line to a factor variable: https://stackoverflow.com/a/16350805/5792088
  facet_wrap(~cohort + parent) +
  theme_classic()
p

####################################################################################
## Question: do older fathers have a higher fraction of early RT variants?
####################################################################################

full_df = bind_rows(ilmn_rt_df, pb_rt_df, decode_rt_df)
full_df %>% group_by(cohort, parent) %>% tally
full_phased_df = full_df %>%  filter(parent != 'Unphased')

age_df_pcgc = read_tsv('longreadclustersequencing/data/parental_age_at_conception_allPCGC_nonNA.txt')
age_df_pcgc %<>%
  mutate(cohort = ifelse(ID %in% pb_rt_df$ID, 'pb_10', 'None')) %>% 
  mutate(cohort = ifelse(ID %in% ilmn_indel_df$ID, 'ilmn_305', cohort)) %>% 
  filter(cohort != 'None') %>% 
  gather(key = 'parent', value = 'parent_age', -ID, -cohort) %>% 
  mutate(parent = ifelse(parent == 'Maternal_age_at_conception', 'Mother', 'Father'))
age_df = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/age_df.txt')
age_df_3gen = age_df %>% filter(cohort == 'decode_3gen') %>% 
  mutate(parent = ifelse(parent == 'mom', 'Mother', 'Father')) %>% 
  rename(parent_age = parental_age)
age_df = bind_rows(age_df_pcgc, age_df_3gen)

age_df %>% group_by(cohort, parent) %>% tally
full_phased_df %>% group_by(cohort) %>% tally
full_phased_df %>% select(ID, cohort) %>% unique %>% group_by(cohort) %>% tally
full_phased_df$Phase %>% unique

cell_phase_group_order = c('Early (G1b/S1)', 'Mid (S2/S3)', 'Late (S4/G2)')
i_tiles = 5
p = full_phased_df %>% 
  group_by(parent, cohort) %>% mutate(parent_age_tile = ntile(parent_age, i_tiles)) %>% ungroup %>% 
  # filter(parent == 'Father') %>% 
  filter(Indel_Class == 'CCC') %>% ## comment for all
  # mutate(cell_phase_group = ifelse(Phase %in% c('G1b', 'S1'), 'Early (G1b/S1)', 'Late (S4/G2)')) %>% 
  # mutate(cell_phase_group = ifelse(Phase %in% c('S2', 'S3'), 'Mid (S2/S3)', cell_phase_group)) %>% 
  # mutate(cell_phase_group = factor(cell_phase_group, levels = cell_phase_group_order)) %>% 
  mutate(Phase = factor(Phase, levels = phase_order)) %>% 
  # filter(!(Phase %in% c('G1b', 'G2'))) %>%
  ## , color = cell_phase_group
  ggplot(aes(x = parent_age_tile, fill = Phase)) +
  geom_histogram(position = 'fill', bins = i_tiles, alpha = 1) +
  # geom_density(alpha = 0.2) +
  facet_wrap(~parent + cohort) +
  scale_fill_manual(values = c('red', 'red3', 'darkred', 'grey50', 'grey70', 'grey90')) +
  ## c('red', 'grey60', 'grey80')
  ylab('Fraction of CCC indels within bin') + xlab('Parental age (quintile)') + #ggtitle('CCC') +
  theme_classic()
  # group_by(cohort, parent, Indel_Class, Phase, parent_age_tile) %>% tally
p
ggsave('longreadclustersequencing/indel_analysis/rep_timing_results/fraction_alltimes_sd1_ccc.png',
       p, width = 6, height = 3.5)

full_phased_df %>% 
  filter(Indel_Class == 'CCC') %>% ## comment for all
  group_by(parent, cohort) %>% mutate(parent_age_tile = ntile(parent_age, i_tiles)) %>% ungroup %>%
  mutate(cell_phase_group = ifelse(Phase %in% c('G1b', 'S1'), 'Early (G1b/S1)', 'Late (S4/G2)')) %>%
  mutate(cell_phase_group = ifelse(Phase %in% c('S2', 'S3'), 'Mid (S2/S3)', cell_phase_group)) %>%
  mutate(cell_phase_group = factor(cell_phase_group, levels = cell_phase_group_order)) %>%
  group_by(cohort, parent, Indel_Class, Phase, parent_age_tile) %>%
  summarise(n_indels = n()) %>% 
  ungroup %>% 
  ## ilmn_305 decode_3gen
  filter(cohort == 'decode_3gen', parent == 'Father') %>% 
  # spread(key = cell_phase_group, value = n_indels) %>% 
  # select(`Early (G1b/S1)`:`Late (S4/G2)`) %>%
  spread(key = Phase, value = n_indels) %>% 
  select(G1b:S4) %>%
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>% as.matrix %>% 
  chisq.test


indels_per_id_per_phase = full_phased_df %>% 
  group_by(ID, parent, parent_age, Indel_Class, Phase) %>% summarise(n_indel_per_phase = n()) %>% 
  ungroup

indels_per_id_per_phase_ALL = full_phased_df %>% 
  group_by(ID, parent, parent_age, Phase) %>% summarise(n_indel_per_phase = n()) %>% 
  ungroup %>% mutate(Indel_Class = 'All')
indels_per_id_per_phase %<>% bind_rows(indels_per_id_per_phase_ALL)
  
expand_df = expand.grid('Phase' = phase_order,
                        'ID' = unique(age_df$ID),
                        'parent' = c('Father', 'Mother'),
                        # cohort = unique(indels_per_phase$cohort),
                        'Indel_Class' = class_order) %>% inner_join(age_df)

indels_per_id_per_phase %<>% select(-parent_age) %>% 
  full_join(expand_df) %>% unique %>% 
  mutate(n_indel_per_phase = ifelse(is.na(n_indel_per_phase), 0, n_indel_per_phase))

indels_per_id_per_phase %<>% 
  group_by(ID, parent, parent_age, Indel_Class) %>% 
  mutate(n_indels_per_id = sum(n_indel_per_phase)) %>% ungroup %>%
  mutate(fraction_indels_per_id = n_indel_per_phase/n_indels_per_id) %>% 
  arrange(ID, parent, parent_age, Indel_Class, Phase)

ilmn_rt_df %>% filter(ID == '1-00392') %>% as.data.frame
indels_per_id_per_phase %>% 
  filter(parent == 'Father') %>% filter(Indel_Class == 'CCC') %>% 
  filter(is.na(fraction_indels_per_id))

indels_per_id_per_phase %>% 
  filter(parent == 'Father') %>% filter(Indel_Class == 'CCC') %>% 
  filter(n_indel_per_phase > 0) %>% group_by(Phase) %>% tally

p = indels_per_id_per_phase %>% 
  filter(parent == 'Father') %>% filter(Indel_Class == 'CCC') %>% 
  ## n_indel_per_phase fraction_indels_per_id
  ggplot(aes(x = parent_age, y = fraction_indels_per_id, color = Phase)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  # geom_errorbar(width = 0.25, position = position_dodge(width=0.95),
  #               color = 'black', alpha = 1, size = 0.25) +
  facet_wrap(~Indel_Class + cohort) +
  theme_classic()
p

indels_per_id_per_phase %>% 
  # filter(!is.na(fraction_indels_per_id)) %>% 
  group_by(cohort, parent, Indel_Class, Phase) %>% 
  summarise(r_frac = cor(fraction_indels_per_id, parent_age),
            p_frac = cor.test(fraction_indels_per_id, parent_age)$p.value,
            r_n = cor(n_indel_per_phase, parent_age),
            p_n = cor.test(n_indel_per_phase, parent_age)$p.value) %>% 
  ungroup %>% 
  arrange(p_n)

# ilmn_rt_df %>% select(contains('ENC')) %>% cor
# ilmn_rt_df %>% 
#   rowwise() %>% 
#   mutate(median_rt = mean(c(ENCFF001GHF, ENCFF001GHH, ENCFF001GHM, ENCFF001GHO,
#                               ENCFF001GHS, ENCFF001GHV))) %>% 
#   ungroup %>% 
#   mutate(pat_age_quint = ntile(Paternal_age_at_conception, 5)) %>% 
#   # mutate(Indel_Class = factor(Indel_Class)) %>% 
#   # aov(median_rt ~ Indel_Class, .) %>% summary
#   ### 1 quantile is youngest
#   group_by(Indel_Class, pat_age_quint) %>% 
#   # summarise(rt_median = median(median_rt), rt_mean = mean(median_rt))
#   summarise_at(vars(contains('ENC')), median)






