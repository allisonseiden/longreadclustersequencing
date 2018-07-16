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

# make dictionary with key: chromosome and value: list of de novos in that chromosome
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
    # initialize dictionary for bounds, key will be de novo and value will be list of bounds for that de novo
    bounds = {};
    # collect correct phased vcf file
    vcf_df = d[chromosome];
    # collect list of de novos for current chromosome
    de_novo_list = dnvs[chromosome];
    # parse through each de novo in current chromosome
    for dn in de_novo_list:
        # initialize bounds list for this de novo variant
        bounds[dn] = [];
        # find index of de novo variant in vcf file
        dnv_index = vcf_df.index[vcf_df['POS'] == dn].item();
        # collect haplotype of de novo variant
        hap = vcf_df['1-00801'][dnv_index];
        # initially set upper discontinuity bound to de novo index to start
        u_discon= dnv_index;
        # move up the file until a discontinuity is found
        while hap[:3] != "0/1" and hap[:3] != "1/0":
            u_discon -= 1;
            hap = vcf_df['1-00801'][u_discon];
        # add upper bound to list for de novo
        bounds[dn].append(vcf_df['POS'][u_discon]);
        # reassign current haplotype to be de novo haplotype
        hap = vcf_df['1-00801'][dnv_index];
        # initially set lower discontinuity bound to de novo index to start
        l_discon = dnv_index;
        # move down the file until a discontinuity is found
        while hap[:3] != "0/1" and hap[:3] != "1/0":
            l_discon += 1;
            hap = vcf_df['1-00801'][l_discon];
        # add lower bound to list for de novo
        bounds[dn].append(vcf_df['POS'][l_discon]);
    return bounds;

# function call to create dictionary with key: chromosome, value: dictionary of de novos and their bounds
all_bounds = {};
for chr in dnvs:
    vcf_df = d[chr];
    all_bounds[chr] = search_discon(vcf_df, chr);


to_phase = {};
for chr in all_bounds:
    vcf_df = d[chr];
    for dnv in all_bounds[chr]:
        # de novo is dnv
        # list of bounds for de novo is all_bounds[chr][dnv]
        curr_bounds = all_bounds[chr][dnv];
        to_phase[dnv] = [];
        u_index = vcf_df.index[vcf_df['POS'] == curr_bounds[0]].item();
        l_index = vcf_df.index[vcf_df['POS'] == curr_bounds[1]].item();
        position = u_index;
        while position <= l_index:
            child = vcf_df['1-00801'][position];
            mom = vcf_df['1-00801-01'][position];
            dad = vcf_df['1-00801-02'][position];
            if child[:3] == "0|1" or child[:3] == "1|0":
                if mom[:3] == "0/0" and (dad[:3] == "1/1" or dad[:3] == "0/1"):
                    to_phase[dnv].append(vcf_df['POS'][position]);
                if dad[:3] == "0/0" and (mom[:3] == "1/1" or mom[:3] == "0/1"):
                    to_phase[dnv].append(vcf_df['POS'][position]);
            position += 1;


phased_to_parent = {};
for dnv in phase:
    phased_to_parent[dnv] = [];
    dnv_index = vcf_df.index[vcf_df['POS'] == dnv].item();
    de_novo_hap = vcf_df['1-00801'][dnv_index];
    for var in phase[dnv]:
        index = vcf_df.index[vcf_df['POS'] == var].item();
        if len(vcf_df['REF'][index]) > 1 or len(vcf_df['ALT'][index]) > 1:
            print("Skipped variant at position " + str(vcf_df['POS'][index]));
            continue;
        child = vcf_df['1-00801'][index];
        mom = vcf_df['1-00801-01'][index];
        dad = vcf_df['1-00801-02'][index];
        if child[:3] == de_novo_hap[:3]:
            if mom[1:3] == "/1":
                phased_to_parent[dnv].append("mom");
            else:
                phased_to_parent[dnv].append("dad");
        if child[:3] != de_novo_hap[:3]:
            if mom[1:3] == "/1":
                phased_to_parent[dnv].append("dad");
            else:
                phased_to_parent[dnv].append("mom");
