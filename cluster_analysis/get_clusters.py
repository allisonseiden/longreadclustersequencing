import pandas as pd
import numpy as np

df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/phasing_analysis/phasing_analysis_df.txt',
                    sep='\t');

clusters = (df['Closest DNV Distance'] != 0) & (df['Closest DNV Distance'] < 20000);

clusters = clusters.astype(int);

df['Cluster'] = clusters;

length = df.shape[0];
for i in range(length):
    if df['Cluster'][i] == 1:
        info = str(df['ID'][i]) + '\t' + str(df['Chrom'][i]) + '\t' + str(df['Location'][i]) + '\t' + str(df['Closest DNV Distance'][i]);
        print(info);
