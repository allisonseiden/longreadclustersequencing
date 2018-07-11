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

dnv_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
parent_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

analysis_df = dnv_df.join(parent_df, how='left');


ti_series = (((analysis_df['Ref'] == 'A') & (analysis_df['Alt'] == 'G')) |
             ((analysis_df['Ref'] == 'G') & (analysis_df['Alt'] == 'A')) |
             ((analysis_df['Ref'] == 'C') & (analysis_df['Alt'] == 'T')) |
             ((analysis_df['Ref'] == 'T') & (analysis_df['Alt'] == 'C')));
tv_series = (((analysis_df['Ref'] == 'A') & ((analysis_df['Alt'] == 'T') | (analysis_df['Alt'] == 'C'))) |
             ((analysis_df['Ref'] == 'G') & ((analysis_df['Alt'] == 'T') | (analysis_df['Alt'] == 'C'))) |
             ((analysis_df['Ref'] == 'T') & ((analysis_df['Alt'] == 'A') | (analysis_df['Alt'] == 'G'))) |
             ((analysis_df['Ref'] == 'C') & ((analysis_df['Alt'] == 'A') | (analysis_df['Alt'] == 'G'))));


analysis_df['Ti'] = ti_series;
analysis_df['Tv'] = tv_series;
analysis_df['Ti'] = analysis_df['Ti'].astype(int);
analysis_df['Tv'] = analysis_df['Tv'].astype(int);

analysis_df = analysis_df[['Ref', 'Alt', 'Ti', 'Tv', 'From Mom', 'From Dad', 'Unphased']];

analysis_df.reset_index(level='Location', inplace=True);

def find_difference(group):
    loc_list = group['Location'].tolist();
    length = len(loc_list);
    distance_list = [];
    for i in range(0, length):
        if length == 1:
            distance_list.append(0);
        elif i == 0:
            d = abs(loc_list[i+1] - loc_list[i]);
            distance_list.append(d);
        elif i == length - 1:
            d = abs(loc_list[i] - loc_list[i-1]);
            distance_list.append(d);
        else:
            d_1 = abs(loc_list[i] - loc_list[i-1]);
            d_2 = abs(loc_list[i+1] - loc_list[i]);
            if (d_1 <= d_2):
                distance_list.append(d_1);
            else:
                distance_list.append(d_2);
    distance_series = pd.Series(data=distance_list);
    return distance_series;



grouped = analysis_df.groupby(['ID', 'Chrom']);
# grouped_loc = grouped['Location'];
print(grouped.apply(find_difference));
