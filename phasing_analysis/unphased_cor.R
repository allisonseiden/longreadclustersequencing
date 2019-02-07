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
# decode data
########################################

decode_indels = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/classified_indels.txt')
decode_dnvs = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/decode_DNMs.tsv')
age_df = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/age_df.txt')


indel_class_vec = c('CCC', 'non-CCC', 'HR')
unphased_ids = decode_dnvs %>% filter(!(Proband_nr %in% age_df_3gen$ID)) %>% 
  select(Proband_nr) %>% unique %>% unlist %>% as.character
age_df_wide = age_df %>% spread(key = parent, value = parental_age)
expand_df = expand.grid('ID' = unphased_ids,
                        # 'ins_del' = c('insertion', 'deletion'),
                        'Indel_Class' = indel_class_vec) %>%
  inner_join(age_df_wide) %>% unique
expand_df %>% group_by(cohort, Indel_Class) %>% tally
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

fit = unphased_indel_ct %>% group_by(cohort, Indel_Class) %>% ## ins_del
  do(model = lm(Indel_ct ~ mom + dad, .))
fit %>% tidy(model) %>% filter(term != '(Intercept)')

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


class_order = c('HR', 'CCC', 'non-CCC', 'All')
p = unphased_indel_ct %>% 
  mutate(Indel_Class = factor(Indel_Class, levels = class_order)) %>% 
  ggplot(aes(x = dad, y = Indel_ct)) +
  geom_point(color = 'grey70', size = 0.25) +
  geom_smooth(method = 'lm', color = 'dodgerblue') +
  facet_wrap(~Indel_Class, nrow = 1) +
  xlab('Paternal age') + ylab('Indels per trio') +
  theme_classic()
p
ggsave('longreadclustersequencing/indel_analysis/unphased_decode_cor.png',
       p, width = 6, height = 2)

dunphased_indel_ct %>% 
  # filter(is.na(dad) | is.na(mom))
  group_by(Indel_Class) %>% 
  summarise(r_dad = cor(Indel_ct, dad),
            r_mom = cor(Indel_ct, mom),
            # r_kendall = cor(Indel_ct, parental_age, method = 'kendall'),
            # r_spearman = cor(Indel_ct, parental_age, method = 'spearman'),
            p_dad = cor.test(Indel_ct, dad, method = 'pearson')$p.value,
            p_mom = cor.test(Indel_ct, mom, method = 'pearson')$p.value,
            # ci_lo = cor.test(Indel_ct, parental_age)$conf.int[[1]],
            # ci_hi = cor.test(Indel_ct, parental_age)$conf.int[[2]],
            Indel_ct = sum(Indel_ct),
            n = n()) %>% ungroup


