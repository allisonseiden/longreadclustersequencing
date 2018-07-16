import pandas as pd
import numpy as np

analysis_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                            sep='\t');
analysis_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
print(analysis_df);
