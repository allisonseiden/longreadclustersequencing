import pandas as pd
import numpy as np

# read in BED file for sample 1-00801
dnv_100801 = pd.read_table('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/1-00801.hg38.dnv.bed',
                            sep='\t', names = ['Chrom', 'Start', 'End', 'Ref', 'Var', 'ID']);

# find number of de novo variants in sample 1-00801
num_dnv = dnv_100801.shape[0];

# create dictionary of chromosome numbers and their corresponding phased VCF files
d = {};
for i in range(1,23):
    num = str(i);
    d["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/1-00801/1-00801_chr' + num + '_phased.vcf',
                                    sep='\t', names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '1-00801', '1-00801-01', '1-00801-02'],
                                    comment = '#');

# make sure not to use same chromosome's vcf multiple times
chrom_list = [];
for row in range(0, num_dnv):
    chr_num = dnv_100801.loc[row, 'Chrom'];
    if chr_num not in chrom_list:
        chrom_list.append(chr_num);

# make dictionary of chromosome with all the de novos in that chromosome
dnvs = {};
for chrom in chrom_list:
    indices = dnv_100801.index[dnv_100801['Chrom'] == chrom].tolist();
    dnv_list = [];
    for index in indices:
        dnv_list.append(dnv_100801['End'][index]);
    dnvs[chrom] = dnv_list;

# no longer need full bed file, all relevant data placed in dnvs dictionary; free memory
del dnv_100801;

# function to search for discontinuities
def search_discon(vcf_df, chromosome):
    bounds = {};
    vcf_df = d[chromosome];
    de_novo_list = dnvs[chromosome];
    for dn in de_novo_list:
        bounds[dn] = [];
        dnv_index = vcf_df.index[vcf_df['POS'] == dn].item();
        hap = vcf_df['1-00801'][dnv_index];
        u_discon= dnv_index;
        while hap[:3] != "0/1" and hap[:3] != "1/0":
            u_discon -= 1;
            hap = vcf_df['1-00801'][u_discon];
        upper_pos = vcf_df['POS'][u_discon];
        bounds[dn].append(upper_pos);
        hap = vcf_df['1-00801'][dnv_index];
        l_discon = dnv_index;
        while hap[:3] != "0/1" and hap[:3] != "1/0":
            l_discon += 1;
            hap = vcf_df['1-00801'][l_discon];
        lower_pos = vcf_df['POS'][l_discon];
        bounds[dn].append(lower_pos);
    return bounds;

# make function call to create dictionary of all bounds
all_bounds = [];
for chr in dnvs:
    vcf_df = d[chr];
    all_bounds.append(search_discon(vcf_df, chr));
    print(all_bounds);
