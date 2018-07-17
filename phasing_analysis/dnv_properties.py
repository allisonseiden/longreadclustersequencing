import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');


group_by_ref_alt = analysis_df.groupby(['ID', 'Ref', 'Alt']);
C_A_df_list = [];


for name, group in group_by_ref_alt:
    if name[1] == 'C' and name[2] == 'A':
        C_A_df_list.append(group);

C_A_df = pd.concat(C_A_df_list, ignore_index=True);
C_A_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

C_A_df = C_A_df[['PB_Mom', 'PB_Dad', 'PB_Unphased', 'IL_Mom', 'IL_Dad', 'IL_Unphased']];

C_A_df = C_A_df.sum(numeric_only=True);
print(C_A_df);
