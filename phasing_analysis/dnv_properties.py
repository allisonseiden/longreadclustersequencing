import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
total_data = analysis_df.sum(numeric_only=True);
total_data.drop(['Location', 'Closest DNV Distance'], inplace=True);
total_data = total_data.append(pd.Series(['All'], index=['Mutational_Class']));
print(total_data);

group_by_ref_alt = analysis_df.groupby(['ID', 'Ref', 'Alt']);
C_A_df_list = [];


for name, group in group_by_ref_alt:
    if (name[1] == 'C' and name[2] == 'A') or (name[1] == 'G' and name[2] == 'T'):
        C_A_df_list.append(group);

C_A_df = pd.concat(C_A_df_list, ignore_index=True);
C_A_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

C_A_data = C_A_df.sum(numeric_only=True);
C_A_data.drop(['Closest DNV Distance'], inplace=True);
C_A_data = C_A_data.append(pd.Series(['C->A'], index=['Mutational_Class']));
print(C_A_data);
