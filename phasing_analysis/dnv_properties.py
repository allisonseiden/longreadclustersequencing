#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Analyze DNV properties.

:Authors: Allison Seiden, Felix Richter
:Date: 2019-01-12
:Copyright: 2018, Allison Seiden
:License: CC BY-SA

python3

"""

import pandas as pd
# import numpy as np

# home_dir = ('/Users/allisonseiden/Documents/longreadclustersequencing/' +
#             'phasing_analysis/')
home_dir = ('/Users/felixrichter/Dropbox/PhD/longreadclustersequencing/' +
            'phasing_analysis/')

# read in dataframe with SNPs and dataframe with indels
analysis_df = pd.read_table(home_dir + 'phasing_analysis_df.txt', sep='\t')
indels_df = pd.read_table(home_dir + 'indels_df.txt', sep='\t')


# cols_interest = ['ID', 'Chrom', 'Location', 'Ref', 'Alt', 'PB_Mom', 'PB_Dad',
#                  'PB_Unphased', 'CpG', 'CpG_Island']
# pacbio_df = analysis_df.loc[:, cols_interest]
pacbio_df = analysis_df

indels_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)
pacbio_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)

indels_df.head()
pacbio_df.head()

df = pacbio_df.join(indels_df, how='left')
df.reset_index(inplace=True)
df.fillna(value=0, inplace=True)

# create column that explicity states whether the de novo came from mom or dad,
# helpful for graph later on
length = df.shape[0]
for i in range(length):
    if df['PB_Mom'][i] == 1:
        df.loc[i, 'Parent'] = 'Mother'
    elif df['PB_Dad'][i] == 1:
        df.loc[i, 'Parent'] = 'Father'
    else:
        if df['PB_Mom_Indel'][i] == 1:
            df.loc[i, 'Parent'] = 'Mother'
            df.loc[i, 'PB_Unphased'] = 0
        elif df['PB_Dad_Indel'][i] == 1:
            df.loc[i, 'Parent'] = 'Father'
            df.loc[i, 'PB_Unphased'] = 0
        else:
            df.loc[i, 'Parent'] = 'None'
    df.loc[i, 'Mom'] = df.loc[i, 'PB_Mom'] + df.loc[i, 'PB_Mom_Indel']
    df.loc[i, 'Dad'] = df.loc[i, 'PB_Dad'] + df.loc[i, 'PB_Dad_Indel']
    df.loc[i, 'Unphased'] = (
        df.loc[i, 'PB_Unphased'] + df.loc[i, 'PB_Unphased_Indel'])


df.head()

df = df[['ID', 'Chrom', 'Location', 'Ref', 'Alt', 'CpG', 'Mom', 'Dad', 'Unphased', 'Parent']]



df.to_csv(path_or_buf='all_pacbio_data.txt', sep='\t')

# get totals in order to get percentages later on
dnv_nums_dict = {}
dnv_nums_dict['Total'] = df.sum(numeric_only=True)
dnv_nums_dict['Total'].drop(['Location'], inplace=True)
dnv_nums_dict['Total'] = dnv_nums_dict['Total'].append(pd.Series(['All'], index=['Mutational_Class']))

group_by_ref_alt = df.groupby(['Ref', 'Alt'])
mutation_groups = {'C_A': [], 'C_T': [], 'C_G': [], 'T_A': [], 'T_C': [],
                   'T_G': [], 'CpG_TpG': [], 'Indels': []}
full_data_dfs = {}

group_list = []
# group together all columns of dataframe that are of the same mutational class
for name, group in group_by_ref_alt:
    if (name[0] == 'C' and name[1] == 'A') or (
            name[0] == 'G' and name[1] == 'T'):
        mutation_groups['C_A'].append(group)
        # group['mut_group'] = 'C_A'
    if (name[0] == 'C' and name[1] == 'T') or (
            name[0] == 'G' and name[1] == 'A'):
        mutation_groups['CpG_TpG'].append(group.loc[lambda df: df.CpG == 1, :])
        mutation_groups['C_T'].append(group.loc[lambda df: df.CpG == 0, :])
    if (name[0] == 'C' and name[1] == 'G') or (
            name[0] == 'G' and name[1] == 'C'):
        mutation_groups['C_G'].append(group)
    if (name[0] == 'T' and name[1] == 'A') or (
            name[0] == 'A' and name[1] == 'T'):
        mutation_groups['T_A'].append(group)
    if (name[0] == 'T' and name[1] == 'C') or (
            name[0] == 'A' and name[1] == 'G'):
        mutation_groups['T_C'].append(group)
    if (name[0] == 'T' and name[1] == 'G') or (
            name[0] == 'A' and name[1] == 'C'):
        mutation_groups['T_G'].append(group)
    if len(name[0]) > 1 or len(name[1]) > 1:
        mutation_groups['Indels'].append(group)



pb_total_mom = float(dnv_nums_dict['Total']['Mom'])
pb_total_dad = float(dnv_nums_dict['Total']['Dad'])
pb_total_unphased = float(dnv_nums_dict['Total']['Unphased'])

# loop through mutational classes and obtain percentages for mom and dad
dnv_percent_dict = {}
for key in mutation_groups:
    full_data_dfs[key] = pd.concat(mutation_groups[key], ignore_index=True)
    dnv_nums_dict[key] = full_data_dfs[key].groupby(['Parent']).sum(numeric_only=True)
    pb_mom_series = dnv_nums_dict[key]['Mom']/pb_total_mom
    pb_dad_series = dnv_nums_dict[key]['Dad']/pb_total_dad
    temp = pd.concat([pb_mom_series, pb_dad_series], axis=1)
    temp.reset_index(inplace=True)
    mutational_class_series = pd.Series([key] * 30)
    dnv_percent_dict[key] = pd.concat([temp, mutational_class_series], axis=1)
    dnv_percent_dict[key].rename(columns={0 : 'Mutational_Class'}, inplace=True)
    dnv_percent_dict[key]['Fraction'] = dnv_percent_dict[key].iloc[:, 1:3].sum(axis=1)
    dnv_percent_dict[key] = dnv_percent_dict[key][['Mutational_Class', 'Parent', 'Fraction']]
    dnv_percent_dict[key].drop([2], inplace=True)
    dnv_percent_dict[key] = dnv_percent_dict[key][dnv_percent_dict[key].Fraction != 0]



dnv_data = pd.concat(dnv_percent_dict)
dnv_data.set_index(['Mutational_Class', 'Parent'], inplace=True)



dnv_data.to_csv(path_or_buf='dnv_parent_percentages.txt', sep='\t')
