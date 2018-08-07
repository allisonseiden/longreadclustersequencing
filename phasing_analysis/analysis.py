import pandas as pd
import numpy as np
import pybedtools

""" Program to analyze phasing information given through Whatshap and parent
    assignment code """


patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389', '1-04460',
                '1-04537', '1-05443', '1-05673', '1-05846'];
bed_list = [];
parent_df_list = [];


# Create lists of pandas dataframes for BED files, Pacbio dataframes, and
# Illumina dataframes for all patient IDs
for ID in patientIDs:
    bed_list.append(pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/' + ID + '_dnv.bed',
                                sep='\t', names = ['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']));
    parent_df_list.append(pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/pacbio_dataframes/' + ID + '_dataframe.txt',
                                sep='\t'));


# Create combined dataframe for all BED files and parentally assigned dataframes
dnv_df = pd.concat(bed_list, ignore_index=True);
parent_df = pd.concat(pb_df_list, ignore_index=True);

dnv_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
parent_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);

# Join together BED dataframe and parent dataframe
temp_one_df = dnv_df.join(parent_df, how='left');


# Create series for which de novos are transitions (ti_series) and which
# de novos are transversions (tv_series)
ti_series = (((temp_one_df['Ref'] == 'A') & (temp_one_df['Alt'] == 'G')) |
             ((temp_one_df['Ref'] == 'G') & (temp_one_df['Alt'] == 'A')) |
             ((temp_one_df['Ref'] == 'C') & (temp_one_df['Alt'] == 'T')) |
             ((temp_one_df['Ref'] == 'T') & (temp_one_df['Alt'] == 'C')));
tv_series = (((temp_one_df['Ref'] == 'A') & ((temp_one_df['Alt'] == 'T') | (temp_one_df['Alt'] == 'C'))) |
             ((temp_one_df['Ref'] == 'G') & ((temp_one_df['Alt'] == 'T') | (temp_one_df['Alt'] == 'C'))) |
             ((temp_one_df['Ref'] == 'T') & ((temp_one_df['Alt'] == 'A') | (temp_one_df['Alt'] == 'G'))) |
             ((temp_one_df['Ref'] == 'C') & ((temp_one_df['Alt'] == 'A') | (temp_one_df['Alt'] == 'G'))));


temp_one_df['Ti'] = ti_series;
temp_one_df['Tv'] = tv_series;
temp_one_df['Ti'] = temp_one_df['Ti'].astype(int);
temp_one_df['Tv'] = temp_one_df['Tv'].astype(int);

temp_one_df = temp_one_df[['Ref', 'Alt', 'Ti', 'Tv', 'Mom', 'Dad', 'Unphased']];

temp_one_df.reset_index(level='Location', inplace=True);

# Function to find the difference between bordering de novos in BED file
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

grouped = temp_one_df.groupby(['ID', 'Chrom']);
temp_one_df = grouped.apply(find_difference);

temp_one_df.set_index(['Location'], append=True, inplace=True);

# Find CpG regions
cpg_bed_list = [];
for ID in patientIDs:
    cpg_bed = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/phasing_analysis/get_fasta_bed/' + ID + '_tri.dnv.bed',
                                        sep=':|-|\t', names=['Chrom', 'Start', 'End', 'Tri_Nucleotide'], engine='python');
    cpg_id = [];
    for elem in cpg_bed['Chrom']:
        cpg_id.append(ID);
    cpg_bed['ID'] = cpg_id;
    cpg_bed_list.append(cpg_bed);

cpg_df = pd.concat(cpg_bed_list, ignore_index=True);
dnv_list = [];
for elem in cpg_df['End']:
    dnv_list.append(elem-1);
dnv_series = pd.Series(dnv_list);
cpg_df['Location'] = dnv_series;
cpg_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
cpg_df['Tri_Nucleotide'] = cpg_df['Tri_Nucleotide'].str.upper();
cpg_df = cpg_df[['Tri_Nucleotide']];
cpg = [];
for elem in cpg_df['Tri_Nucleotide']:
    if (elem[1] == 'C' and (elem[0] == 'G' or elem[2] == 'G')) or (elem[1] == 'G' and (elem[0] == 'C' or elem[2] == 'C')):
        cpg.append(1);
    else:
        cpg.append(0);
cpg_df['CpG'] = cpg;

temp_two_df = temp_one_df.join(cpg_df, how='left');

# Find CpG islands using bedtools intersect with bed files and CpGisland file
# from UCSC Genome Browser
dnv_bed_list = [];
for ID in patientIDs:
    dnv_bed = pybedtools.BedTool('/hpc/users/seidea02/longreadclustersequencing/data/' + ID + '_dnv.bed');
    dnv_bed.intersect('CpG_islands.bed').saveas('CpG_islands/CpG_islands_' + ID + '.bed');
    dnv_bed_list.append(pd.read_table('CpG_islands/CpG_islands_' + ID + '.bed', sep='\t',
                        names=['Chrom', 'Start', 'Location', 'Ref', 'Alt', 'ID']));

dnv_bed_df = pd.concat(dnv_bed_list, ignore_index=True);
dnv_bed_df.set_index(['ID', 'Chrom', 'Location'], inplace=True);
dnv_bed_df['CpG_Island'] = [True] * dnv_bed_df.shape[0];
dnv_bed_df = dnv_bed_df[['CpG_Island']];
print(dnv_bed_df);

analysis_df = temp_two_df.join(dnv_bed_df, how='left');
analysis_df.fillna(value=False, inplace=True);
analysis_df['CpG_Island'] = analysis_df['CpG_Island'].astype(int);

analysis_df.to_csv(path_or_buf='phasing_analysis_df.txt', sep='\t')
print(analysis_df);
