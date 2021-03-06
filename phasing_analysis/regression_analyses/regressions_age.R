# Felix Richter
# 1/22/2019
# Indel correlations with parental age
#################################################

setwd('D:/Dropbox/PhD/')
setwd('/Users/frichter/Dropbox (Personal)/PhD/')
setwd('/Users/felixrichter/Dropbox/PhD/')
options(stringsAsFactors=FALSE)

p = c('magrittr', 'extraDistr', 'broom', 'car', 'lme4',
      'purrr', 'dplyr', 'ggplot2', 'tidyr', 'readr', 'wesanderson', 'stringi')
lapply(p, require, character.only = TRUE)

# prepared in regression_prep.R
indel_full_ct = read_tsv('longreadclustersequencing/indel_analysis/counts_per_id_ilmn_pb_repeat_HR_AT_GC_2019_02_06.txt')
# prepared in decode_indel_correlations.R
indel_decode = read_tsv('longreadclustersequencing/literature/jonsson_decode_2017/indel_cts_per_id_repeat_HRATGC_2019_02_06.txt')

# counts_per_id_ilmn_pb_repeat_HR_AT_GC_2019_02_06.txt
# counts_per_id_ilmn_pb_2019_01_22.txt
# indel_cts_per_id_repeat_HRATGC_2019_02_06.txt
# indel_cts_per_id.txt
####################################
# Join with fraction phased
####################################

# prepared in  fraction_phased.R
fract_phased = read_tsv('longreadclustersequencing/phasing_analysis/fraction_phased/fraction_phased_2019_01_22.txt')

# now only keep the fraction of SNVs phased
snv_phased = fract_phased %>% filter(indel != 'indel') %>% select(-indel) %>% 
  rename(fraction_snv_phased = fraction_phased)
indel_full_ct %<>% inner_join(snv_phased)

# confirm all important columns have values
indel_full_ct %>% bind_rows(indel_decode) %>% 
  filter(is.na(parental_age), parental_age == 0, is.na(Indel_ct),
         is.na(fraction_snv_phased))

# join your data with DECODE
indel_full_ct %<>% bind_rows(indel_decode) 

# confirm all indel classes have every cohort ID accounted for
indel_full_ct %>% group_by(cohort, parent, Indel_Class) %>% tally %>% as.data.frame

# indel_full_ct %>% write_tsv('longreadclustersequencing/indel_analysis/full_cts_all_2019_03_03.txt')

####################################
# Run correlations and plot
####################################

# Sum over repeats for HR AT vs GC
cor_tbl = indel_full_ct %>% 
  filter(grepl('HR', Indel_Class)) %>%
  group_by(cohort, Indel_Class, parent, ID, parental_age) %>% 
  summarise(Indel_ct = sum(Indel_ct)) %>% ungroup %>% 
  group_by(cohort, Indel_Class, parent) %>% 
  summarise(r = cor(Indel_ct, parental_age),
            # r_kendall = cor(Indel_ct, parental_age, method = 'kendall'),
            r_spearman = cor(Indel_ct, parental_age, method = 'spearman'),
            p = cor.test(Indel_ct, parental_age, method = 'spearman')$p.value,
            # ci_lo = cor.test(Indel_ct, parental_age)$conf.int[[1]],
            # ci_hi = cor.test(Indel_ct, parental_age)$conf.int[[2]],
            Indel_ct = sum(Indel_ct),
            n = n()) %>% ungroup

# cor_tbl %>% write_tsv('longreadclustersequencing/phasing_analysis/regression_analyses/cor_tbl_HR_AT_GC.txt')

# check the correlations
cor_tbl = indel_full_ct %>% 
  group_by(cohort, Indel_Class, parent, in_repeat) %>% ## in_repeat
  summarise(r = cor(Indel_ct, parental_age, method = 'spearman'), ## pearson
            # r_kendall = cor(Indel_ct, parental_age, method = 'kendall'),
            # r_spearman = cor(Indel_ct, parental_age, method = 'spearman'),
            p = cor.test(Indel_ct, parental_age, method = 'spearman', exact = F)$p.value,
            delta = 1.96*(1/(sqrt(n() - 3))), 
            ci_lo = tanh(atanh(r) - delta),
            ci_hi = tanh(atanh(r) + delta),
            # ci_lo = cor.test(Indel_ct, parental_age)$conf.int[[1]],
            # ci_hi = cor.test(Indel_ct, parental_age)$conf.int[[2]],
            Indel_ct = sum(Indel_ct),
            n = n()) %>% ungroup %>% select(-delta)

cor_tbl
# cor_tbl %>% 
#   mutate_at(vars(r, p, ci_lo, ci_hi), function(x) signif(x, digits = 3)) %>% 
#   write_tsv('longreadclustersequencing/phasing_analysis/regression_analyses/cor_tbl.txt')
# cor_tbl_repeat.txt

####################################
# calculate meta-analysis p-values
####################################

cor_sub = cor_tbl %>% filter(Indel_Class == 'CCC', parent == 'dad') ## in_repeat == 'N', 
m1 = meta::metacor(cor_sub$r, cor_sub$n, sm='ZCOR')
summary(m1)
meta::forest(m1)

cor_tbl %>% #filter(grepl('decode', cohort)) %>% 
  mutate(transformed_p = -2*log(p)) %>% filter(!is.na(transformed_p)) %>% 
  ## transform the correlation coefficients
  # mutate(transformed_r_to_z = 0.5*log((1+r)/(1-r))) %>% 
  # mutate(study_r_se = 1/(sqrt(n - 3))) %>% head %>% as.data.frame
  arrange(cohort, in_repeat) %>% 
  # filter(Indel_Class == 'CCC') %>% 
  group_by(Indel_Class, parent, in_repeat) %>% ## in_repeat
  summarise(df_i = 2*n(), p_meta = 1 - pchisq(sum(transformed_p), df = df_i),
            p_list = paste(signif(p, 2), collapse = ', '),
            # cohort_list = paste(cohort, collapse = ', '),
            indel_sum = sum(Indel_ct),
            Indel_ct = paste(Indel_ct, collapse = ', '), 
            r_list = paste(signif(r, 2), collapse = ', '))

####################################
# plotting correlations
####################################

class_order = c('HR', 'CCC', 'non-CCC', 'All') %>% rev
cohort_order = c('deCODE (N=225)\nHaplotype-sharing',
                 'Illumina (N=305)\nShort-read tracing',
                 'PacBio (N=10)\nLong-read tracing')
plot_y_min = -0.1

p = cor_tbl %>% 
  mutate(ci_lo = ifelse(ci_lo < plot_y_min, plot_y_min, ci_lo)) %>% 
  mutate(Indel_Class = factor(Indel_Class, levels = class_order)) %>% 
  ## clean the cohort names (potentially)
  mutate(cohort = ifelse(cohort == 'decode_3gen', cohort_order[[1]], cohort)) %>%
  mutate(cohort = ifelse(cohort == 'ilmn_305', cohort_order[[2]], cohort)) %>%
  mutate(cohort = ifelse(cohort == 'pb_10', cohort_order[[3]], cohort)) %>%
  mutate(cohort = factor(cohort, levels = rev(cohort_order))) %>% 
  ggplot(., aes(x = Indel_Class, y = r, ymin = ci_lo, ymax = ci_hi, fill = parent, color = parent)) +
  geom_col(width = 0.85, position = 'dodge2', alpha = 0.75) +
  scale_fill_manual(values = c('#352846', '#E2AE38')) +
  scale_color_manual(values = c('#352846', '#E2AE38')) +
  geom_errorbar(width = 0.25, position = position_dodge(width=0.85),
                color = 'black', alpha = 1, size = 0.25) +
  facet_wrap(~cohort) + 
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(plot_y_min, 1)) +
  coord_flip() +
  theme_classic()
p
ggsave('longreadclustersequencing/phasing_analysis/regression_analyses/cor_bar_plots_long.png',
       p, width = 6, height = 1.95)

####################################
# plotting for repeats
####################################


# class_order = c('HR', 'CCC', 'non-CCC', 'All') %>% rev
cohort_order = c('deCODE (N=225)\nHaplotype-sharing',
                 'Illumina (N=305)\nShort-read tracing',
                 'PacBio (N=10)\nLong-read tracing')
plot_y_min = -0.1

p = cor_tbl %>% 
  filter(Indel_Class == 'CCC') %>% 
  mutate(ci_lo = ifelse(ci_lo < plot_y_min, plot_y_min, ci_lo)) %>%
  mutate(in_repeat = ifelse(in_repeat == 'Y', 'In repeat', 'Not repeat')) %>% 
  ## clean the cohort names (potentially)
  mutate(cohort = ifelse(cohort == 'decode_3gen', cohort_order[[1]], cohort)) %>%
  mutate(cohort = ifelse(cohort == 'ilmn_305', cohort_order[[2]], cohort)) %>%
  mutate(cohort = ifelse(cohort == 'pb_10', cohort_order[[3]], cohort)) %>%
  mutate(cohort = factor(cohort, levels = rev(cohort_order))) %>% 
  ggplot(., aes(x = in_repeat, y = r, ymin = ci_lo, ymax = ci_hi, fill = parent, color = parent)) +
  geom_col(width = 0.95, position = 'dodge2', alpha = 0.75) +
  scale_fill_manual(values = c('#352846', '#E2AE38')) +
  scale_color_manual(values = c('#352846', '#E2AE38')) +
  geom_errorbar(width = 0.25, position = position_dodge(width=0.95),
                color = 'black', alpha = 1, size = 0.25) +
  facet_wrap(~cohort) + 
  xlab('') + ylab('') +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(plot_y_min, 1)) +
  coord_flip() +
  theme_classic()
p
# ggsave('longreadclustersequencing/phasing_analysis/regression_analyses/cor_bar_plots_repeat.png',
#        p, width = 6, height = 1.55)

####################################
# Plot scatterplots
####################################
cohort_order = c('deCODE (N=225)\nHaplotype-sharing',
                 'Illumina (N=305)\nShort-read tracing',
                 'PacBio (N=10)\nLong-read tracing')
p = indel_full_ct %>% 
  # mutate(ci_lo = ifelse(ci_lo < plot_y_min, plot_y_min, ci_lo)) %>% 
  # mutate(Indel_Class = factor(Indel_Class, levels = class_order)) %>% 
  filter(Indel_Class == 'All') %>% 
  ## clean the cohort names (potentially)
  mutate(cohort = ifelse(cohort == 'decode_3gen', cohort_order[[1]], cohort)) %>%
  mutate(cohort = ifelse(cohort == 'ilmn_305', cohort_order[[2]], cohort)) %>%
  mutate(cohort = ifelse(cohort == 'pb_10', cohort_order[[3]], cohort)) %>%
  mutate(cohort = factor(cohort, levels = rev(cohort_order))) %>%
  ggplot(., aes(x = parental_age, y = Indel_ct, fill = parent, color = parent)) +
  geom_point(alpha = 0.5, size = 0.75) +
  geom_smooth(method = 'lm', se = F, alpha = 0.25) +
  scale_fill_manual(values = c('#352846', '#E2AE38')) +
  scale_color_manual(values = c('#352846', '#E2AE38')) +
  facet_wrap(~cohort) + 
  xlab('') + ylab('') +
  theme_classic()
p
ggsave('longreadclustersequencing/phasing_analysis/regression_analyses/indel_cor_all.png',
       p, width = 6.5, height = 2.25)
# indel_cor_nonccc.png indel_cor_ccc.png indel_cor_hr.png

####################################################
# looking at the relative burden
####################################################

dad_age = indel_full_ct %>% filter(parent == 'dad') %>% 
  select(cohort, ID, parental_age) %>% unique
fraction_paternal = indel_full_ct %>% 
  select(-parental_age, -fraction_snv_phased) %>% 
  spread(key = parent, value = Indel_ct) %>% 
  # group_by(cohort, Indel_Class) %>% tally ## confirm sample sizes are correct
  mutate(fraction_paternal = dad/(dad + mom)) %>% 
  inner_join(dad_age)

fraction_paternal %>% 
  filter(Indel_Class == 'All') %>% 
  ggplot(., aes(x = parental_age, y = fraction_paternal)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~cohort+Indel_Class) +
  theme_classic()

fraction_paternal %>% 
  group_by(Indel_Class, cohort) %>% 
  summarise(r = cor(fraction_paternal, parental_age, use = 'complete.obs'),
            r_kendall = cor(fraction_paternal, parental_age, method = 'kendall', use = 'complete.obs'),
            r_spearman = cor(fraction_paternal, parental_age, method = 'spearman', use = 'complete.obs'),
            p = cor.test(fraction_paternal, parental_age, method = 'pearson', use = 'complete.obs')$p.value,
            # ci_lo = cor.test(Indel_ct, parental_age)$conf.int[[1]],
            # ci_hi = cor.test(Indel_ct, parental_age)$conf.int[[2]],
            n = n()) %>% ungroup
  # filter(Indel_Class == 'All') %>% 
  

####################################################
# re-run correlations as linear regressions
####################################################

fit_df = indel_full_ct %>% 
  group_by(cohort, Indel_Class, parent, in_repeat) %>% 
  ## confirm that standardizing gives you the R^2 values. It does!
  mutate(Indel_ct_std = (Indel_ct - mean(Indel_ct))/sd(Indel_ct),
         parental_age_std = (parental_age - mean(parental_age))/sd(parental_age),
         fraction_snv_phased_std = (fraction_snv_phased - mean(fraction_snv_phased))/sd(fraction_snv_phased)) %>% 
  do(model = lm(Indel_ct ~ parental_age, .)) ## fraction_snv_phased
  # do(model = lm(Indel_ct_std ~ parental_age_std, .)) ## + fraction_snv_phased_std
  # do(model = glm(Indel_ct ~ parental_age, family = 'poisson', .))

fit_df %>% tidy(model, conf.int = T) %>% 
  # write_tsv('longreadclustersequencing/phasing_analysis/fraction_phased/multiple_regression.txt')
  filter(term != '(Intercept)') %>% 
  # filter(cohort == 'pb_10')
  # filter(term == 'fraction_snv_phased') %>%
  arrange(p.value)
fit_df %>% glance(model)

# include the fraction phased and ancestry as covariates in regression

ilmn_indel_full_ct %>% head

################################################
# Trying to bootstrap the CI since it's >= 1
# https://stats.stackexchange.com/a/175029/227447
################################################
