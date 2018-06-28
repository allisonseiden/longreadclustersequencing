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

# for chrom in chrom_list:
#     print(chrom)
#     hap = d[chrom]['1-00801'];
#     print(hap);



# def search_discon(self, chromosome):
#     hap = d[chromosome]['1-00801'];


# #print(chr22_100801.loc[:,['POS','1-00801']]);
# pos_hap = chr22[['POS', '1-00801']];
# hap = chr22['1-00801'];
# length = hap.size;


# for i in range(0,length):
#     mod_hap = hap[i][:3];
#     if mod_hap == "0/1" or mod_hap == "1/0":
#         print(mod_hap);
