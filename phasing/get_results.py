import pandas as pd
import numpy as np

full_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt', sep='\t', engine='python');

illumina_df = full_df.loc[:, ['ID', 'Chrom', 'Location', 'IL_Mom', 'IL_Dad', 'IL_Unphased']]


illumina_df['IL_Phased'] = illumina_df['IL_Mom'] + illumina_df['IL_Dad']
illumina_df = illumina_df[['ID', 'IL_Phased', 'IL_Unphased']];

illumina_df = illumina_df.groupby(['ID']).sum()
illumina_df.reset_index(inplace=True);

illumina_df['Fraction'] = illumina_df['IL_Phased']/(illumina_df['IL_Phased'] + illumina_df['IL_Unphased']);

il_std_fractions = illumina_df['Fraction'].std();
il_perc_phased = illumina_df['IL_Phased'].sum()/(illumina_df['IL_Phased'].sum() + illumina_df['IL_Unphased'].sum());

print(il_perc_phased);
print(il_std_fractions)


pacbio_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/phasing_analysis/all_pacbio_data.txt', sep='\t', engine='python')

pacbio_df['Phased'] = pacbio_df['Mom'] + pacbio_df['Dad']
pacbio_df = pacbio_df[['ID', 'Phased', 'Unphased']];

pacbio_df = pacbio_df.groupby(['ID']).sum()
pacbio_df.reset_index(inplace=True);

pacbio_df['Fraction'] = pacbio_df['Phased']/(pacbio_df['Phased'] + pacbio_df['Unphased']);

pb_std_fractions = pacbio_df['Fraction'].std();
pb_perc_phased = pacbio_df['Phased'].sum()/(pacbio_df['Phased'].sum() + pacbio_df['Unphased'].sum());

print(pb_perc_phased);
print(pb_std_fractions)
