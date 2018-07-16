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


temp_one_df = dnv_df.join(pb_parent_df, how='left');
temp_two_df = temp_one_df.join(il_parent_df, how='left');


ti_series = (((temp_two_df['Ref'] == 'A') & (temp_two_df['Alt'] == 'G')) |
             ((temp_two_df['Ref'] == 'G') & (temp_two_df['Alt'] == 'A')) |
             ((temp_two_df['Ref'] == 'C') & (temp_two_df['Alt'] == 'T')) |
             ((temp_two_df['Ref'] == 'T') & (temp_two_df['Alt'] == 'C')));
tv_series = (((temp_two_df['Ref'] == 'A') & ((temp_two_df['Alt'] == 'T') | (temp_two_df['Alt'] == 'C'))) |
             ((temp_two_df['Ref'] == 'G') & ((temp_two_df['Alt'] == 'T') | (temp_two_df['Alt'] == 'C'))) |
             ((temp_two_df['Ref'] == 'T') & ((temp_two_df['Alt'] == 'A') | (temp_two_df['Alt'] == 'G'))) |
             ((temp_two_df['Ref'] == 'C') & ((temp_two_df['Alt'] == 'A') | (temp_two_df['Alt'] == 'G'))));


temp_two_df['Ti'] = ti_series;
temp_two_df['Tv'] = tv_series;
temp_two_df['Ti'] = temp_two_df['Ti'].astype(int);
temp_two_df['Tv'] = temp_two_df['Tv'].astype(int);

temp_two_df = temp_two_df[['Ref', 'Alt', 'Ti', 'Tv', 'PB_Mom', 'PB_Dad', 'PB_Unphased', 'IL_Mom', 'IL_Dad', 'IL_Unphased']];

temp_two_df.reset_index(level='Location', inplace=True);

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

grouped = temp_two_df.groupby(['ID', 'Chrom']);
temp_two_df = grouped.apply(find_difference);

temp_two_df.set_index(['Location'], append=True, inplace=True);

cpg_bed_list = [];
for ID in patientIDs:
    cpg_bed_list.append(pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/get_fasta_bed/' + ID + '_tri.hg38.dnv.bed',
                                        sep=':|-|\t', names=['Chrom', 'Start', 'End', 'Tri_Nucleotide', 'ID'], engine='python'));

cpg_df = pd.concat(cpg_bed_list, ignore_index=True);
print(cpg_df);




dnv_bed_list = [];
for ID in patientIDs:
    dnv_bed = pybedtools.BedTool('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/' + ID + '.hg38.dnv.bed');
    dnv_bed.intersect('CpG_islands.bed').saveas('CpG_islands/CpG_islands_' + ID + '.bed');
    dnv_bed_list.append(pd.read_table('CpG_islands/CpG_islands_' + ID + '.bed', sep='\t',
                        names=['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']));

dnv_bed_df = pd.concat(dnv_bed_list, ignore_index=True);
CpG_i_length = dnv_bed_df.shape[0];
all_ones = [1] * CpG_i_length;
CpG_island = pd.Series(all_ones, dtype=int);
print(CpG_island);
dnv_bed_df['CpG_Island'] = CpG_island;
dnv_bed_df = dnv_bed_df[['ID', 'Chrom', 'Location', 'CpG_Island']];
dnv_bed_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

analysis_df = temp_two_df.join(dnv_bed_df, how='left');
analysis_df.fillna(value=0, inplace=True);



# print(analysis_df);
