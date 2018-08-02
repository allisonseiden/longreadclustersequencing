import pandas as pd
import numpy as np

indels_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/phasing_analysis/indels_df.txt',
                            sep='\t');
classified_indels_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/indel_analysis/classified_indels.txt',
                            sep='\t', names=['ID', 'Chrom', 'Start', 'Location', 'Ref',
                                                'Alt', 'Allele', 'Indel_Class',
                                                'repName', 'repClass', 'repFamily']);

classified_indels_df = classified_indels_df[['ID', 'Chrom', 'Location', 'Ref',
                                              'Alt', 'Allele', 'Indel_Class',
                                              'repName', 'repClass', 'repFamily']];

indels_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
classified_indels_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

df = indels_df.join(classified_indels_df, how='left');
df.reset_index(inplace=True);
df.to_csv(path_or_buf='all_indel_info.txt', sep='\t');

length = df.shape[0];
for i in range(length):
    if df['PB_Mom_Indel'][i] == 1:
        df.loc[i, 'Parent'] = 'Mother';
    elif df['PB_Dad_Indel'][i] == 1:
        df.loc[i, 'Parent'] = 'Father';
    else:
        df.loc[i, 'Parent'] = 'None';

dnv_nums_dict = {};
dnv_nums_dict['Total'] = df.sum(numeric_only=True);
dnv_nums_dict['Total'].drop(['Location'], inplace=True);
dnv_nums_dict['Total'] = dnv_nums_dict['Total'].append(pd.Series(['All'], index=['Indel_Class']));


group_by_class = df.groupby(['Indel_Class']);
indel_classes = {'HR' : [], 'CCC' : [], 'non-CCC' : []};
full_data_dfs = {};

for name, group in group_by_class:
    if name == 'HR':
        indel_classes['HR'].append(group);
    if name == 'CCC':
        indel_classes['CCC'].append(group);
    if name == 'non-CCC':
        indel_classes['non-CCC'].append(group);


total_mom = float(dnv_nums_dict['Total']['PB_Mom_Indel']);
total_dad = float(dnv_nums_dict['Total']['PB_Dad_Indel']);
total_unphased = float(dnv_nums_dict['Total']['PB_Unphased_Indel']);

dnv_percent_dict = {};
for key in indel_classes:
    full_data_dfs[key] = pd.concat(indel_classes[key], ignore_index=True);
    dnv_nums_dict[key] = full_data_dfs[key].groupby(['Parent']).sum(numeric_only=True);
    mom_series = dnv_nums_dict[key]['PB_Mom_Indel']/total_mom;
    dad_series = dnv_nums_dict[key]['PB_Dad_Indel']/total_dad;
    temp = pd.concat([mom_series, dad_series], axis=1);
    temp.reset_index(inplace=True);
    indel_class_series = pd.Series([key] * 30);
    dnv_percent_dict[key] = pd.concat([temp, indel_class_series], axis=1);
    dnv_percent_dict[key].rename(columns={0 : 'Indel_Class'}, inplace=True);
    dnv_percent_dict[key]['Fraction'] = dnv_percent_dict[key].iloc[:, 1:3].sum(axis=1);
    dnv_percent_dict[key] = dnv_percent_dict[key][['Indel_Class', 'Parent', 'Fraction']];
    dnv_percent_dict[key].drop([2], inplace=True);
    dnv_percent_dict[key] = dnv_percent_dict[key][dnv_percent_dict[key].Fraction != 0];


dnv_data = pd.concat(dnv_percent_dict);
dnv_data.set_index(['Indel_Class', 'Parent'], inplace=True);



dnv_data.to_csv(path_or_buf='dnv_parent_percentages.txt', sep='\t');
