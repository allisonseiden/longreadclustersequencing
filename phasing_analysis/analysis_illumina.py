"""Analyze phasing information from Whatshap and parent assignment code.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-12-06
:Copyright: 2018, Allison Seiden, Felix Richter
:License: CC BY-SA

cd /hpc/users/richtf01/longreadclustersequencing/data/gmkf2

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix/0.2.6 bedtools/2.27.1
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
python3


"""

import pandas as pd
import glob
import re
# import numpy as np
# import pybedtools


home_dir = '/hpc/users/richtf01/longreadclustersequencing/'

"""Prepare DNV dataframe."""

dnv_iter = glob.iglob(home_dir + 'data/gmkf2/1-*_dnv.bed')
dnv_col_names = ['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']
dnv_list = []
for dnv_f in dnv_iter:
    dnv_list.append(pd.read_table(dnv_f, sep='\t', names=dnv_col_names))


dnv_df = pd.concat(dnv_list, ignore_index=True)
dnv_df.head()
dnv_df.shape

"""Prepare assigned DNV dataframe."""

assigned_dir = home_dir + 'phasing_analysis/illumina_dataframes_2018_11_25/'
assigned_iter = glob.iglob(assigned_dir + '1-*_dataframe.txt')
assigned_list = []
for assigned_f in assigned_iter:
    assigned_list.append(pd.read_table(assigned_f, sep='\t'))


assigned_df = pd.concat(assigned_list, ignore_index=True)
assigned_df.head()
assigned_df.shape

"""Combine full DNV call info with parental assignment info."""

dnv_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)
assigned_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)

# Join together BED dataframe and parent dataframe
dnv_full_df = dnv_df.join(assigned_df, how='left')


"""Transition and transversion info."""
# Create series for which de novos are transitions (ti_series) and which
# de novos are transversions (tv_series)
ti_series = (((dnv_full_df['Ref'] == 'A') & (dnv_full_df['Alt'] == 'G')) |
             ((dnv_full_df['Ref'] == 'G') & (dnv_full_df['Alt'] == 'A')) |
             ((dnv_full_df['Ref'] == 'C') & (dnv_full_df['Alt'] == 'T')) |
             ((dnv_full_df['Ref'] == 'T') & (dnv_full_df['Alt'] == 'C')))
tv_series = (((dnv_full_df['Ref'] == 'A') & (
                (dnv_full_df['Alt'] == 'T') | (dnv_full_df['Alt'] == 'C'))) |
             ((dnv_full_df['Ref'] == 'G') & (
                (dnv_full_df['Alt'] == 'T') | (dnv_full_df['Alt'] == 'C'))) |
             ((dnv_full_df['Ref'] == 'T') & (
                (dnv_full_df['Alt'] == 'A') | (dnv_full_df['Alt'] == 'G'))) |
             ((dnv_full_df['Ref'] == 'C') & (
                (dnv_full_df['Alt'] == 'A') | (dnv_full_df['Alt'] == 'G'))))


dnv_full_df['Ti'] = ti_series
dnv_full_df['Tv'] = tv_series
dnv_full_df['Ti'] = dnv_full_df['Ti'].astype(int)
dnv_full_df['Tv'] = dnv_full_df['Tv'].astype(int)


"""Keep columns of interest and find distance to closest DNV."""
dnv_full_df = dnv_full_df[['Ref', 'Alt', 'Ti', 'Tv', 'Mom', 'Dad', 'Unphased']]
dnv_full_df.reset_index(level='Location', inplace=True)


def find_difference(group):
    """Find the difference between neighboring DNVs in same ID."""
    loc_list = group['Location'].tolist()
    length = len(loc_list)
    distance_list = []
    for i in range(0, length):
        if length == 1:
            distance_list.append(0)
        elif i == 0:
            d = abs(loc_list[i+1] - loc_list[i])
            distance_list.append(d)
        elif i == length - 1:
            d = abs(loc_list[i] - loc_list[i-1])
            distance_list.append(d)
        else:
            d_1 = abs(loc_list[i] - loc_list[i-1])
            d_2 = abs(loc_list[i+1] - loc_list[i])
            if (d_1 <= d_2):
                distance_list.append(d_1)
            else:
                distance_list.append(d_2)
    d_series = pd.Series(data=distance_list)
    group['Closest DNV Distance'] = d_series.values
    return group


grouped = dnv_full_df.groupby(['ID', 'Chrom'])
dnv_full_df = grouped.apply(find_difference)

dnv_full_df.set_index(['Location'], append=True, inplace=True)
dnv_full_df.head()
dnv_full_df.shape

"""CpG from trinucleotide context"""

cpg_bed_list = []
cpg_iter = glob.iglob(home_dir + 'data/gmkf2/trinuc_out/1-*.txt')
cpgcol = ['Chrom', 'Start', 'End', 'Tri_Nucleotide']
for cpg_f in cpg_iter:
    # ID = re.sub('.*CpG_islands/|_dnv_cpg_.*', '', cpg_f)
    ID = re.sub('.*trinuc_out/|_dnv_trinuc.*', '', cpg_f)
    cpg_bed = pd.read_table(cpg_f, sep=':|-|\t', names=cpgcol, engine='python')
    # add a column for ID
    cpg_bed['ID'] = [ID for i in cpg_bed['Chrom']]
    cpg_bed_list.append(cpg_bed)

cpg_df = pd.concat(cpg_bed_list, ignore_index=True)
# subtract 1 from the end position, add as new column
cpg_df['Location'] = pd.Series([i-1 for i in cpg_df['End']])
cpg_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)
# convert trinucleotide content to uppercase
cpg_df['Tri_Nucleotide'] = cpg_df['Tri_Nucleotide'].str.upper()
cpg_df = cpg_df[['Tri_Nucleotide']]
cpg = []
for elem in cpg_df['Tri_Nucleotide']:
    if (elem[1] == 'C' and (elem[0] == 'G' or elem[2] == 'G')) or (
            elem[1] == 'G' and (elem[0] == 'C' or elem[2] == 'C')):
        cpg.append(1)
    else:
        cpg.append(0)


cpg_df['CpG'] = cpg

# if happy, change temp_df to dnv_full_df
dnv_full_df = dnv_full_df.join(cpg_df, how='left')
dnv_full_df.head()
dnv_full_df.shape

"""Cpg island overlap"""
cpg_isle_bed_list = []
cpg_isle_iter = glob.iglob(home_dir + 'data/gmkf2/CpG_islands/1-*.bed')
cpg_cols = ['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']
for cpg_isle_f in cpg_isle_iter:
    cpg_isle_bed = pd.read_table(cpg_isle_f, sep='\t', names=cpg_cols)
    cpg_isle_bed_list.append(cpg_isle_bed)

cpg_isle_bed_df = pd.concat(cpg_isle_bed_list, ignore_index=True)
cpg_isle_bed_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)
cpg_isle_bed_df['CpG_Island'] = [True] * cpg_isle_bed_df.shape[0]
cpg_isle_bed_df = cpg_isle_bed_df[['CpG_Island']]
print(cpg_isle_bed_df)

analysis_df = dnv_full_df.join(cpg_isle_bed_df, how='left')
analysis_df.fillna(value=False, inplace=True)
analysis_df['CpG_Island'] = analysis_df['CpG_Island'].astype(int)

analysis_df.head()
analysis_df.shape

out_f = home_dir + 'phasing_analysis/phasing_analysis_df_ilmn_2018_12_06.txt'
analysis_df.to_csv(path_or_buf=out_f, sep='\t')


# save dataframe here
"""Subset just the indels and join those"""


#
