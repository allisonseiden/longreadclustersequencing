import pandas as pd
import numpy as np

analysis_df = pd.read_table('/Users/Seiden/Documents/Summer2018/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
dnv_data_dict = {};
dnv_data_dict['Total'] = analysis_df.sum(numeric_only=True);
dnv_data_dict['Total'].drop(['Location', 'Ti', 'Tv', 'Closest DNV Distance'], inplace=True);
dnv_data_dict['Total'] = dnv_data_dict['Total'].append(pd.Series(['All'], index=['Mutational_Class']));

group_by_ref_alt = analysis_df.groupby(['Ref', 'Alt']);
mutation_groups = {'C_A' : [], 'C_T' : [], 'C_G' : [], 'T_A' : [], 'T_C' : [], 'T_G' : [], 'CpG_TpG' : []};
# mutation_groups = {'C_A' : [], 'G_T' : [], 'C_T' : [], 'G_A' : [], 'C_G' : [], 'G_C' : [],
#                     'T_A' : [], 'A_T' : [], 'T_C' : [], 'A_G' : [], 'T_G' : [], 'A_C' : []};
full_data_dfs = {};

for name, group in group_by_ref_alt:
    if (name[0] == 'C' and name[1] == 'A') or (name[0] == 'G' and name[1] == 'T'):
        mutation_groups['C_A'].append(group);
    if (name[0] == 'C' and name[1] == 'T') or (name[0] == 'G' and name[1] == 'A'):
        mutation_groups['CpG_TpG'].append(group.loc[lambda df: df.CpG == 1, :]);
        mutation_groups['C_T'].append(group.loc[lambda df: df.CpG == 0, :]);
    if (name[0] == 'C' and name[1] == 'G') or (name[0] == 'G' and name[1] == 'C'):
        mutation_groups['C_G'].append(group);
    if (name[0] == 'T' and name[1] == 'A') or (name[0] == 'A' and name[1] == 'T'):
        mutation_groups['T_A'].append(group);
    if (name[0] == 'T' and name[1] == 'C') or (name[0] == 'A' and name[1] == 'G'):
        mutation_groups['T_C'].append(group);
    if (name[0] == 'T' and name[1] == 'G') or (name[0] == 'A' and name[1] == 'C'):
        mutation_groups['T_G'].append(group);
    # if (name[0] == 'C' and name[1] == 'A'):
    #     mutation_groups['C_A'].append(group);
    # if (name[0] == 'G' and name[1] == 'T'):
    #     mutation_groups['G_T'].append(group);
    # if (name[0] == 'C' and name[1] == 'T'):
    #     mutation_groups['C_T'].append(group);
    # if (name[0] == 'G' and name[1] == 'A'):
    #     mutation_groups['G_A'].append(group);
    # if (name[0] == 'C' and name[1] == 'G'):
    #     mutation_groups['C_G'].append(group);
    # if (name[0] == 'G' and name[1] == 'C'):
    #     mutation_groups['G_C'].append(group);
    # if (name[0] == 'T' and name[1] == 'A'):
    #     mutation_groups['T_A'].append(group);
    # if (name[0] == 'A' and name[1] == 'T'):
    #     mutation_groups['A_T'].append(group);
    # if (name[0] == 'T' and name[1] == 'C'):
    #     mutation_groups['T_C'].append(group);
    # if (name[0] == 'A' and name[1] == 'G'):
    #     mutation_groups['A_G'].append(group);
    # if (name[0] == 'T' and name[1] == 'G'):
    #     mutation_groups['T_G'].append(group);
    # if (name[0] == 'A' and name[1] == 'C'):
    #     mutation_groups['A_C'].append(group);

pb_total_mom = float(dnv_data_dict['Total']['PB_Mom']);
pb_total_dad = float(dnv_data_dict['Total']['PB_Dad']);
pb_total_phased = float(dnv_data_dict['Total']['PB_Mom'] + dnv_data_dict['Total']['PB_Dad']);
pb_total_unphased = float(dnv_data_dict['Total']['PB_Unphased']);
il_total_mom = float(dnv_data_dict['Total']['IL_Mom']);
il_total_dad = float(dnv_data_dict['Total']['IL_Dad']);
il_total_phased = float(dnv_data_dict['Total']['IL_Mom'] + dnv_data_dict['Total']['IL_Dad']);
il_total_unphased = float(dnv_data_dict['Total']['IL_Unphased']);

for key in mutation_groups:
    full_data_dfs[key] = pd.concat(mutation_groups[key], ignore_index=True);
    full_data_dfs[key].set_index(['ID', 'Chrom', 'Location'], inplace=True);
    dnv_data_dict[key] = full_data_dfs[key].sum(numeric_only=True);
    dnv_data_dict[key].drop(['Ti', 'Tv', 'Closest DNV Distance'], inplace=True);
    pb_phased = (dnv_data_dict[key]['PB_Mom'] + dnv_data_dict[key]['PB_Dad'])/pb_total_phased;
    pb_mom = dnv_data_dict[key]['PB_Mom']/pb_total_mom;
    pb_dad = dnv_data_dict[key]['PB_Dad']/pb_total_dad;
    pb_unphased = (dnv_data_dict[key]['PB_Unphased'])/pb_total_unphased;
    il_phased = (dnv_data_dict[key]['IL_Mom'] + dnv_data_dict[key]['IL_Dad'])/il_total_phased;
    il_mom = dnv_data_dict[key]['IL_Mom']/il_total_mom;
    il_dad = dnv_data_dict[key]['IL_Dad']/il_total_dad;
    il_unphased = (dnv_data_dict[key]['IL_Unphased'])/il_total_unphased;
    dnv_data_dict[key] = dnv_data_dict[key].append(pd.Series([pb_mom, pb_dad, pb_phased, pb_unphased, il_mom, il_dad, il_phased, il_unphased], index=['PB_Percent_Mom', 'PB_Percent_Dad', 'PB_Percent_Phased', 'PB_Percent_Unphased', 'IL_Percent_Mom', 'IL_Percent_Dad', 'IL_Percent_Phased', 'IL_Percent_Unphased']));
    dnv_data_dict[key] = dnv_data_dict[key].append(pd.Series([str(key)], index=['Mutational_Class']));

dnv_data = pd.DataFrame(dnv_data_dict).transpose();
dnv_data.reset_index(inplace=True);
dnv_data.set_index(['Mutational_Class'], inplace=True);
dnv_data.drop(['index'], axis=1, inplace=True);

dnv_data = dnv_data[['PB_Mom', 'PB_Percent_Mom', 'PB_Dad', 'PB_Percent_Dad', 'PB_Percent_Phased',
                    'PB_Unphased', 'PB_Percent_Unphased', 'IL_Mom', 'IL_Percent_Mom',
                    'IL_Dad', 'IL_Percent_Dad', 'IL_Percent_Phased', 'IL_Unphased',
                    'IL_Percent_Unphased', 'CpG', 'CpG_Island']];


dnv_data.to_csv(path_or_buf='dnv_data.txt', sep='\t');
