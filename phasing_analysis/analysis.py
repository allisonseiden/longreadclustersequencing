import pandas as pd
import numpy as np

# last patient ID not in list because fixing bam
patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389', '1-04460',
                '1-04537', '1-05443', '1-05673'];
bed_list = [];
df_list = [];


for ID in patientIDs:
    bed_list.append(pd.read_table('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/' + ID + '.hg38.dnv.bed',
                                sep='\t', names = ['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']));
    df_list.append(pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/' + ID + '_dataframe.txt',
                                sep='\t'));

dnv_df = pd.concat(bed_list, ignore_index=True);
parent_df = pd.concat(df_list, ignore_index=True);

print(parent_df);

# analysis_df = dnv_df.join(parent_df, on=['ID', 'Chrom', 'Location'], how='left');
# print(analysis_df);
