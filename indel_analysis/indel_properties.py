import pandas as pd
import numpy as np

patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389', '1-04460', '1-04537', '1-05443', '1-05673', '1-05846'];

indels_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/phasing_analysis/indels_df.txt',
                            sep='\t');
classified_indels_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/indel_analysis/classified_indels.txt',
                            sep='\t')

classified_indels_df.rename(columns={'End' : 'Location'}, inplace=True)

classified_indels_df = classified_indels_df[['ID', 'Chrom', 'Location', 'Ref',
                                              'Alt', 'Allele', 'HR', 'CCC', 'non-CCC',
                                              'repName', 'repClass', 'repFamily']];


indels_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
classified_indels_df.sort_values(by=['ID', 'Chrom', 'Location'], inplace=True);
classified_indels_df.set_index(['ID', 'Chrom', 'Location'], inplace=True)


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

# ----------- get totals instead of fractions

# df = df[['ID', 'PB_Mom_Indel', 'PB_Dad_Indel', 'PB_Unphased_Indel', 'HR', 'CCC', 'non-CCC', 'Parent']];
# group_list = [];
# group_by_ID = df.groupby(['ID'])
# for name, group in group_by_ID:
#     parent_list = group['Parent'].tolist();
#     if 'Mother' not in parent_list:
#         group = group.append({'ID' : name, 'Parent' : 'Mother', 'PB_Mom_Indel' : 0, 'PB_Dad_Indel' : 0, 'PB_Unphased_Indel' : 0, 'HR' : 0, 'CCC' : 0, 'non-CCC' : 0}, ignore_index=True);
#     if 'Father' not in parent_list:
#         group = group.append({'ID' : name, 'Parent' : 'Father', 'PB_Mom_Indel' : 0, 'PB_Dad_Indel' : 0, 'PB_Unphased_Indel' : 0, 'HR' : 0, 'CCC' : 0, 'non-CCC' : 0}, ignore_index=True);
#     group_list.append(group);
#
# dnv_nums = pd.concat(group_list, ignore_index=True);
# dnv_sums = dnv_nums.groupby(['ID', 'Parent']).sum();
# dnv_sums['Num_DNVs'] = dnv_sums.iloc[:, 0:2].sum(axis=1);
# dnv_sums.reset_index(inplace=True);
# dnv_sums = dnv_sums[dnv_sums['Parent'] != 'None']
# dnv_sums = dnv_sums[['ID', 'Parent', 'Num_DNVs', 'HR', 'CCC', 'non-CCC']]
# dnv_sums.to_csv(path_or_buf='dnv_parent_totals_by_ID.txt', sep='\t');

# dnv_nums_dict['Total'] = df.sum(numeric_only=True);
# dnv_nums_dict['Total'].drop(['Location'], inplace=True);
# dnv_nums_dict['Total'] = dnv_nums_dict['Total'].append(pd.Series(['All'], index=['Indel_Class']));

# ------- get fractions

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
    dnv_nums_dict[key] = full_data_dfs[key].groupby(['ID']).sum(numeric_only=True);
    mom_series = dnv_nums_dict[key]['PB_Mom_Indel']/total_mom;
    dad_series = dnv_nums_dict[key]['PB_Dad_Indel']/total_dad;
    temp = pd.concat([mom_series, dad_series], axis=1);
    temp.reset_index(inplace=True);
    indel_class_series = pd.Series([key] * 10);
    dnv_percent_dict[key] = pd.concat([temp, indel_class_series], axis=1);
    dnv_percent_dict[key].rename(columns={0 : 'Indel_Class'}, inplace=True);
    dnv_percent_dict[key]['Fraction'] = dnv_percent_dict[key].iloc[:, 1:3].sum(axis=1);
    dnv_percent_dict[key] = dnv_percent_dict[key][['ID', 'Indel_Class', 'Parent', 'Fraction']];
    dnv_percent_dict[key].drop([2], inplace=True);
    dnv_percent_dict[key] = dnv_percent_dict[key][dnv_percent_dict[key].Num_DNVs != 0];


dnv_data = pd.concat(dnv_percent_dict);
dnv_data.set_index(['ID', 'Indel_Class', 'Parent'], inplace=True);



dnv_data.to_csv(path_or_buf='dnv_parent_percentages_by_ID.txt', sep='\t');
