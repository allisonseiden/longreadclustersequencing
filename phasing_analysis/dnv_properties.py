import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
analysis_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

group_by_ref_alt = analysis_df.groupby(['ID', 'Ref', 'Alt']);


for name, group in group_by_ref_alt:
    if name[1] == 'C' and name[2] == 'A':
        print(type(group));
