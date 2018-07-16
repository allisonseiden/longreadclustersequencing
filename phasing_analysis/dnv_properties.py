import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
phased_unphased_pacbio = analysis_df.groupby(['ID', 'Chrom', 'Location',
                                                        'PB_Mom', 'PB_Dad',
                                                        'PB_Unphased']);
for name, group in phased_unphased_pacbio:
    print(name);
    print(group);               
