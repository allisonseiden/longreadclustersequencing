import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
analysis_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

C_A_pb_phased = ((analysis_df['Ref'] == 'C') & (analysis_df['Alt'] == 'A') &
              ((analysis_df['PB_Mom'] == 1) | (analysis_df['PB_Dad'] == 1)));

C_A_il_phased = ((analysis_df['Ref'] == 'C') & (analysis_df['Alt'] == 'A') &
              ((analysis_df['IL_Mom'] == 1) | (analysis_df['IL_Dad'] == 1)));

C_A_phased_unphased = pd.concat([C_A_pb_phased, C_A_il_phased]);

print(C_A_phased_unphased);
