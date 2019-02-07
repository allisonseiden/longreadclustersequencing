# Felix Richter
# 10/22/2018
# Make separate DNV BED files for every trio
######################################################


write_id = function(ID, dn_for_bed) {
  out_f = paste0('longreadclustersequencing/data/gmkf2/gmkf2_w_cause/', ID, '_dnv.bed')
  print(out_f)
  dn_for_bed %>% 
    filter(Blinded.ID %in% ID) %>% 
    write_tsv(out_f, col_names = F)
}

chr_vec = c(paste0('chr', 1:22), 'chrX', 'chrY')

dn_qual %>% group_by(cases_known_cause) %>% tally

dn_for_bed = dn_qual %>% 
  filter(wave == "GMKF2") %>% 
  mutate(Chrom = factor(Chrom, levels = chr_vec)) %>% arrange(Chrom, End) %>% 
  select(Chrom, Start, End, Ref, Alt, Blinded.ID) %>% mutate(Start = Start - 1) %>% 
  mutate_all(as.character) 

per_id_list = map(unique(dn_for_bed$Blinded.ID), write_id, dn_for_bed)

########################################
# Prepare for sorting_hat input
########################################
# dn_qual form enrichment_global.R
chr_vec = c(paste0('chr', 1:22), 'chrX', 'chrY')
dn_bed_df = dn_qual %>% 
  filter(snv_indel == 'indel') %>% 
  mutate(Chrom = factor(Chrom, levels = chr_vec)) %>% arrange(Chrom, End) %>% 
  select(Chrom, Start, End, Ref, Alt, Blinded.ID) %>% mutate(Start = Start - 1) %>% 
  mutate_all(as.character) %>% unique

write_tsv(dn_bed_df, 'longreadclustersequencing/data/dnvs_2019_02_07.bed', col_names = F)


