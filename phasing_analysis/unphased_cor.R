# Felix Richter
# 11/18/2018
# Unphased indel correlations with parental age
#################################################

setwd('D:/Dropbox/PhD/')
setwd('/Users/frichter/Dropbox (Personal)/PhD/')
setwd('/Users/felixrichter/Dropbox/PhD/')
options(stringsAsFactors=FALSE)

p = c('car', 'magrittr', 'extraDistr', 'ggplot2', 'tidyr', 'broom', 
      'readr', 'wesanderson', 'purrr', 'dplyr', 'stringi')
lapply(p, require, character.only = TRUE)

## preparing all DNVs for sorting-hat input in get_dnv_dfs.R

########################################
# SFARI and PCGC data
########################################

dnvs_sorted = read_tsv('longreadclustersequencing/data/dnvs_2019_02_07_dnv_sorting_hat.txt')

## get previously phased IDs (keep only unphased IDs)
phased_ids = read_tsv('longreadclustersequencing/data/phased_ids_pcgc.txt') %>% 
  unlist %>% as.character

dnvs_sorted %<>% filter(!(ID %in% phased_ids))

## get parental age
age_df_cases = read_tsv('longreadclustersequencing/data/parental_age_at_conception_allPCGC_nonNA.txt')

# dn = readRDS(file = "whole_genome/pcgcV2_trios/dn_allTSS_UpAndDownTSS_718mCHD_RefGeneUCSC_19_01_27.RDS")
# age_df_ctrls = dn %>% filter(case_ctrl == 'Ctrl') %>%
#   select(Blinded.ID, Maternal.Age.at.Proband.Birth, Paternal.Age.at.Proband.Birth) %>% unique
# names(age_df_ctrls) = names(age_df_cases)
# age_df_ctrls %>% filter(is.na(Maternal_age_at_conception) | is.na(Paternal_age_at_conception))
## account for gestation
# age_df_ctrls %>% 
#   mutate(Maternal_age_at_conception = Maternal_age_at_conception - 0.75,
#          Paternal_age_at_conception = Paternal_age_at_conception - 0.75) %>% 
#   write_tsv('longreadclustersequencing/data/parental_age_at_conception_allSSC_nonNA.txt')

age_df_ctrls = read_tsv('longreadclustersequencing/data/parental_age_at_conception_allSSC_nonNA.txt')
age_df_ilmn = bind_rows('pcgc' = age_df_cases, 'ssc' = age_df_ctrls, .id = 'cohort')

no_par_age_ids = dnvs_sorted %>% anti_join(age_df_ilmn) %>% select(ID) %>% unique %>% unlist %>% as.character

dnvs_sorted %<>% inner_join(age_df_ilmn) 

# dnvs_sorted %>% 
#   mutate(Indel_Class = factor(Indel_Class, levels = c('HR', 'CCC', 'non-CCC'))) %>% 
#   mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% unique %>% 
#   group_by(cohort, Indel_Class, ins_del) %>% tally %>% ungroup %>% 
#   arrange(cohort, Indel_Class, desc(ins_del)) %>% 
#   write_tsv('longreadclustersequencing/indel_analysis/unphased_summary_cts_pcgc.txt')

### 
indel_class_vec = c('CCC', 'non-CCC', 'HR')
unphased_ids = age_df_ilmn %>% filter(!(ID %in% phased_ids)) %>% select(ID) %>% unlist %>% as.character
expand_df = expand.grid('ID' = unphased_ids,
                        # 'ins_del' = c('insertion', 'deletion'),
                        'Indel_Class' = indel_class_vec) %>%
  inner_join(age_df_ilmn) %>% unique
expand_df %>% group_by(cohort, Indel_Class) %>% tally

## actually count indels/trio
unphased_indel_ct = dnvs_sorted %>% 
  mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% unique %>% 
  group_by(ID, Indel_Class) %>% summarise(Indel_ct = n()) %>% ungroup %>% 
  right_join(expand_df) %>% 
  mutate(Indel_ct = ifelse(is.na(Indel_ct), 0, Indel_ct)) 

## include 'All'
unphased_indel_all_df = unphased_indel_ct %>%
  group_by(ID, Paternal_age_at_conception, Maternal_age_at_conception, cohort) %>%
  summarise(Indel_ct = sum(Indel_ct)) %>% ungroup %>% 
  mutate(Indel_Class = 'All')
unphased_indel_ct %<>% bind_rows(unphased_indel_all_df)

unphased_indel_ct %>% group_by(cohort, Indel_Class) %>% tally
unphased_indel_ct_pcgc = unphased_indel_ct  ## needed for combining w decode
fit = unphased_indel_ct %>%
  group_by(Indel_Class, cohort) %>% ## ins_del cohort
  # do(model = lm(Indel_ct ~ Paternal_age_at_conception + Maternal_age_at_conception, .))
  do(model = glm(Indel_ct ~ Paternal_age_at_conception + Maternal_age_at_conception, ., family = 'poisson'))
  # do(model = lme4::glmer(Indel_ct ~ Paternal_age_at_conception + 
  #                          (Paternal_age_at_conception | cohort) +
  #                          Maternal_age_at_conception +
  #                          (Maternal_age_at_conception | cohort), data = .,
  #                        family = 'poisson'))

fit %>% tidy(model) %>% 
  # filter(term != '(Intercept)')
  filter(grepl('age_at_conception$', term))
fit_pcgc = fit

class_order = c('HR', 'CCC', 'non-CCC', 'All')
p = unphased_indel_ct_pcgc %>% 
  filter(cohort == 'pcgc') %>% 
  mutate(Indel_Class = factor(Indel_Class, levels = class_order)) %>% 
  ggplot(aes(x = Paternal_age_at_conception, y = Indel_ct)) +
  geom_point(color = 'grey70', size = 0.25) +
  # geom_smooth(method = 'lm', color = 'dodgerblue') +
  geom_smooth(method = 'glm', color = 'dodgerblue', method.args = list(family = "poisson")) +
  facet_wrap(~Indel_Class, nrow = 1) +
  xlab('Paternal age') + ylab('Indels per trio') +
  theme_classic()
p
ggsave('longreadclustersequencing/indel_analysis/unphased_cors/unphased_pcgc_poisson.png',
       p, width = 6, height = 2)

########################################
# decode data
########################################

decode_indels = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/classified_indels.txt')
decode_dnvs = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/decode_DNMs.tsv')
age_df = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/age_df.txt')
age_df_3gen = age_df %>% filter(cohort == 'decode_3gen')


indel_class_vec = c('CCC', 'non-CCC', 'HR')
unphased_ids = decode_dnvs %>% filter(!(Proband_nr %in% age_df_3gen$ID)) %>% 
  select(Proband_nr) %>% unique %>% unlist %>% as.character
age_df_wide = age_df %>% spread(key = parent, value = parental_age)
expand_df = expand.grid('ID' = unphased_ids,
                        # 'ins_del' = c('insertion', 'deletion'),
                        'Indel_Class' = indel_class_vec) %>%
  inner_join(age_df_wide) %>% unique
expand_df %>% group_by(cohort, Indel_Class) %>% tally ## 1323 
expand_df %>% filter(is.na(dad) | is.na(mom))
unphased_indel_ct = decode_indels %>% 
  mutate(Indel_Class = ifelse(grepl('HR', Indel_Class), 'HR', Indel_Class)) %>% 
  mutate(ins_del = ifelse((stri_length(Ref) > 1), 'deletion', 'insertion')) %>% unique %>% 
  group_by(ID, Indel_Class) %>% summarise(Indel_ct = n()) %>% ungroup %>% 
  right_join(expand_df) %>% 
  mutate(Indel_ct = ifelse(is.na(Indel_ct), 0, Indel_ct)) 

## add all indels
unphased_indel_all_df = unphased_indel_ct %>%
  group_by(ID, dad, mom, cohort) %>%
  summarise(Indel_ct = sum(Indel_ct)) %>% ungroup %>% 
  mutate(Indel_Class = 'All')
unphased_indel_ct %<>% bind_rows(unphased_indel_all_df)

fit = unphased_indel_ct %>%
  group_by(cohort, Indel_Class) %>% ## ins_del
  # do(model = lm(Indel_ct ~ mom + dad, .))
  do(model = glm(Indel_ct ~ mom + dad, ., family = 'poisson'))
fit %>% tidy(model) %>% filter(term != '(Intercept)')
fit_decode = fit

## checking residuals
fit_dad = unphased_indel_ct %>% group_by(cohort, Indel_Class) %>%
  do(model = lm(Indel_ct ~ dad, .))
fit_dad %>% tidy(model) %>% filter(term != '(Intercept)')
fit_dad_res = fit_dad %>% augment(model) %>% as.data.frame %>%
  rename(res_dad = .resid) %>% select(res_dad)
## add in dad's residuals
fit_mom = unphased_indel_ct %>% 
  mutate(res_dad = fit_dad_res$res_dad) %>% 
  group_by(cohort, Indel_Class) %>%
  do(model = lm(res_dad ~ mom, .))
fit_mom %>% tidy(model) %>% filter(term != '(Intercept)')
###

class_order = c('HR', 'CCC', 'non-CCC', 'All')
p = unphased_indel_ct %>% 
  mutate(Indel_Class = factor(Indel_Class, levels = class_order)) %>% 
  ggplot(aes(x = dad, y = Indel_ct)) +
  geom_point(color = 'grey70', size = 0.25) +
  geom_smooth(method = 'glm', color = 'dodgerblue', method.args = list(family = "poisson")) +
  facet_wrap(~Indel_Class, nrow = 1) +
  xlab('Paternal age') + ylab('Indels per trio') +
  theme_classic()
p
ggsave('longreadclustersequencing/indel_analysis/unphased_cors/unphased_decode_poisson.png',
       p, width = 6, height = 2)

unphased_indel_ct_decode = unphased_indel_ct

##########################################
# Combined Poisson regressions
##########################################

##########################################
# Meta-analysis
##########################################

fit_meta = bind_rows(fit_decode %>% tidy(model), fit_pcgc %>% tidy(model)) %>% 
  filter(term != '(Intercept)') %>% 
  rename(parent = term, p_interest = p.value) %>% 
  select(cohort:parent, p_interest, estimate) %>% 
  mutate(parent = ifelse(grepl('Maternal', parent), 'mom', parent)) %>% 
  mutate(parent = ifelse(grepl('Paternal', parent), 'dad', parent))

fit_meta %>% 
  mutate(transformed_p = -2*log(p_interest)) %>% filter(!is.na(transformed_p)) %>% 
  arrange(cohort) %>% 
  group_by(Indel_Class, parent) %>% ## in_repeat
  summarise(df_i = 2*n(), p_meta = 1 - pchisq(sum(transformed_p), df = df_i),
            p_list = paste(signif(p_interest, 2), collapse = ', '),
            # cohort_list = paste(cohort, collapse = ', '),
            # indel_sum = sum(Indel_ct),
            # Indel_ct = paste(Indel_ct, collapse = ', '), 
            beta_list = paste(signif(estimate, 2), collapse = ', '))


##########################################
# Single correlations
# DONT USE THESE BECAUSE DON'T ACCOUNT
# FOR MATERNAL/PATERNAL correlation
##########################################
# 
# cor_tbl_pcgc = unphased_indel_ct_pcgc %>% 
#   group_by(Indel_Class, cohort) %>% 
#   summarise(r_dad = cor(Indel_ct, Paternal_age_at_conception, method = 'spearman'),
#             r_mom = cor(Indel_ct, Maternal_age_at_conception, method = 'spearman'),
#             # r_kendall = cor(Indel_ct, parental_age, method = 'kendall'),
#             # r_spearman = cor(Indel_ct, parental_age, method = 'spearman'),
#             p_dad = cor.test(Indel_ct, Paternal_age_at_conception, method = 'spearman', exact = F)$p.value,
#             p_mom = cor.test(Indel_ct, Maternal_age_at_conception, method = 'spearman', exact = F)$p.value,
#             # ci_lo = cor.test(Indel_ct, parental_age)$conf.int[[1]],
#             # ci_hi = cor.test(Indel_ct, parental_age)$conf.int[[2]],
#             Indel_ct = sum(Indel_ct),
#             n = n()) %>% ungroup
# 
# cor_tbl_decode = unphased_indel_ct_decode %>% 
#   # filter(is.na(dad) | is.na(mom))
#   group_by(Indel_Class, cohort) %>% 
#   summarise(r_dad = cor(Indel_ct, dad, method = 'spearman'),
#             r_mom = cor(Indel_ct, mom, method = 'spearman'),
#             # r_kendall = cor(Indel_ct, parental_age, method = 'kendall'),
#             # r_spearman = cor(Indel_ct, parental_age, method = 'spearman'),
#             p_dad = cor.test(Indel_ct, dad, method = 'spearman', exact = F)$p.value,
#             p_mom = cor.test(Indel_ct, mom, method = 'spearman', exact = F)$p.value,
#             # ci_lo = cor.test(Indel_ct, parental_age)$conf.int[[1]],
#             # ci_hi = cor.test(Indel_ct, parental_age)$conf.int[[2]],
#             Indel_ct = sum(Indel_ct),
#             n = n()) %>% ungroup
# 
# 
# cor_tbl_meta = bind_rows(cor_tbl_decode, cor_tbl_pcgc) 
# 
# cor_tbl_meta %>%
#   mutate(p_interest = p_mom) %>%
#   mutate(transformed_p = -2*log(p_interest)) %>% filter(!is.na(transformed_p)) %>% 
#   arrange(cohort) %>% 
#   group_by(Indel_Class) %>% ## in_repeat
#   summarise(df_i = 2*n(),
#             p_chisq_log = pchisq(sum(transformed_p), df = df_i, log.p = T),
#             p_meta = 1 - pchisq(sum(transformed_p), df = df_i),
#             p_list = paste(signif(p_interest, 2), collapse = ', '),
#             # cohort_list = paste(cohort, collapse = ', '),
#             indel_sum = sum(Indel_ct),
#             Indel_ct = paste(Indel_ct, collapse = ', '),
#             r_list = paste(signif(r_dad, 2), collapse = ', '))


