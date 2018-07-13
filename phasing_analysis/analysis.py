import pandas as pd
import numpy as np
import pybedtools

# last patient ID not in list because fixing bam
patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389', '1-04460',
                '1-04537', '1-05443', '1-05673', '1-05846'];
bed_list = [];
pb_df_list = [];
il_df_list = [];


for ID in patientIDs:
    bed_list.append(pd.read_table('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/' + ID + '.hg38.dnv.bed',
                                sep='\t', names = ['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']));
    pb_df_list.append(pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/pacbio_dataframes/' + ID + '_dataframe.txt',
                                sep='\t'));
    il_df_list.append(pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/illumina_dataframes/' + ID + '_dataframe.txt',
                                sep='\t'));

dnv_df = pd.concat(bed_list, ignore_index=True);
pb_parent_df = pd.concat(pb_df_list, ignore_index=True);
pb_parent_df.rename(columns={"From Mom": "PB_Mom", "From Dad": "PB_Dad", "Unphased": "PB_Unphased"}, inplace=True);
il_parent_df = pd.concat(il_df_list, ignore_index=True);
il_parent_df.rename(columns={"From Mom": "IL_Mom", "From Dad": "IL_Dad", "Unphased": "IL_Unphased"}, inplace=True);

dnv_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
pb_parent_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
il_parent_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

temp_df = dnv_df.join(pb_parent_df, how='left');
analysis_df = temp_df.join(il_parent_df, how='left');



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

analysis_df = analysis_df[['Ref', 'Alt', 'Ti', 'Tv', 'PB_Mom', 'PB_Dad', 'PB_Unphased', 'IL_Mom', 'IL_Dad', 'IL_Unphased']];

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
    d_series = pd.Series(data=distance_list);
    group['Closest DNV Distance'] = d_series.values;
    return group;

grouped = analysis_df.groupby(['ID', 'Chrom']);
analysis_df = grouped.apply(find_difference);

dnv_bed_list = [];
for ID in patientIDs:
    dnv_bed = pybedtools.BedTool('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/' + ID + '.hg38.dnv.bed');
    dnv_bed.intersect('CpG_islands.bed').saveas('CpG_islands/CpG_islands_' + ID + '.bed');
    dnv_bed_list.append(pd.read_table('CpG_islands/CpG_islands_' + ID + '.bed', sep='\t',
                        names=['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']));

dnv_bed_df = pd.concat(dnv_bed_list, ignore_index=True);
dnv_bed_df.set_index(['ID', 'Chrom', 'Location']);
CpG_i_length = dnv_bed_df.shape[0];
all_ones = [1] * CpG_i_length;
CpG_island = pd.Series(all_ones);
dnv_bed_df['CpG_island'] = CpG_island;
analysis_df = analysis_df.join(dnv_bed_df, how='left');
print(analysis_df);
