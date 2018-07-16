import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
analysis_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

# C_A_pb_phased = (((analysis_df['Ref'] == 'C') & (analysis_df['Alt'] == 'A')) &
#               ((analysis_df['PB_Mom']) | (analysis_df['PB_Dad'])));
# C_A_pb_phased = C_A_pb_phased.astype(int);
#
# C_A_il_phased = (((analysis_df['Ref'] == 'C') & (analysis_df['Alt'] == 'A')) &
#               ((analysis_df['IL_Mom']) | (analysis_df['IL_Dad'])));
# C_A_il_phased = C_A_il_phased.astype(int);
group_by_ref_alt = analysis_df.groupby(['Ref', 'Alt']);

# C_A_phased_unphased = pd.concat([C_A_pb_phased, C_A_il_phased], axis=1);

for name, group in group_by_ref_alt:
    print(name);
    print(group);
