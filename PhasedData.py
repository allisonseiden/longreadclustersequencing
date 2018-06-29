import pandas as pd;
import numpy as np;


class PhasedData:
    def __init__(self, patientID):
        self.id = patientID;
        self.mom = patientID + '-01';
        self.dad = patientID + '-02';
        self.bed = pd.read_table('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/' + self.id + '.hg38.dnv.bed',
                                    sep='\t', names = ['Chrom', 'Start', 'End', 'Ref', 'Var', 'ID']);
        self.vcf_dfs = {};
        self.dnvs = {};
        self.bounds = {};
        self.to_phase = {};
        self.phased_to_parent = {};

    def create_vcf_dictionary(self):
        for i in range(1,23):
            num = str(i);
            self.vcf_dfs["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '/' + self.id + '_chr' + num + '_phased.vcf',
                                            sep='\t', names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '1-00801', '1-00801-01', '1-00801-02'],
                                            comment = '#');

    def create_dnvs_dictionary(self):
        num_dnvs = self.bed.shape[0];
        chrom_list = [];
        for row in range(0, num_dnvs):
            chr_num = self.bed.loc[row, 'Chrom'];
            if chr_num not in chrom_list:
                chrom_list.append(chr_num);

        for chrom in chrom_list:
            indices = self.bed.index[self.bed['Chrom'] == chrom].tolist();
            self.dnvs[chrom] = [];
            for index in indices:
                self.dnvs[chrom].append(self.bed['End'][index]);

    def search_discon(self, chromosome):
        chr_bounds = {};
        for chr in self.dnvs:
            curr_vcf = self.vcf_dfs[chr];
            for dnv in self.dnvs[chr]:
                chr_bounds[dnv] = [];
                dnv_index = curr_vcf.index[curr_vcf['POS'] == dnv].item();
                hap = curr_vcf[self.id][dnv_index];
                u_discon= dnv_index;
                while hap[:3] != "0/1" and hap[:3] != "1/0":
                    u_discon -= 1;
                    hap = curr_vcf[self.id][u_discon];
                chr_bounds[dnv].append(curr_vcf['POS'][u_discon]);
                hap = curr_vcf[self.id][dnv_index];
                l_discon = dnv_index;
                while hap[:3] != "0/1" and hap[:3] != "1/0":
                    l_discon += 1;
                    hap = curr_vcf[self.id][l_discon];
                chr_bounds[dnv].append(curr_vcf['POS'][l_discon]);
        return chr_bounds;

    def fill_bounds_dictionary(self):
        for chr in self.dnvs:
            self.bounds[chr] = self.search_discon(chr);

    def find_variants_for_phasing(self):
        for chr in self.bounds:
            curr_vcf = self.vcf_dfs[chr];
            for dnv in self.bounds[chr]:
                # de novo is dnv
                # list of bounds for de novo is all_bounds[chr][dnv]
                curr_bounds = self.bounds[chr][dnv];
                self.to_phase[dnv] = [];
                u_index = int(curr_vcf.index[curr_vcf['POS'] == curr_bounds[0]]);
                l_index = curr_vcf.index[curr_vcf['POS'] == curr_bounds[1]].item();
                position = u_index;
                while position <= l_index:
                    child = curr_vcf[self.id][position];
                    mom = curr_vcf[self.mom][position];
                    dad = curr_vcf[self.dad][position];
                    if child[:3] == "0|1" or child[:3] == "1|0":
                        if mom[:3] == "0/0" and (dad[:3] == "1/1" or dad[:3] == "0/1"):
                            self.to_phase[dnv].append(curr_vcf['POS'][position]);
                        if dad[:3] == "0/0" and (mom[:3] == "1/1" or mom[:3] == "0/1"):
                            self.to_phase[dnv].append(curr_vcf['POS'][position]);
                    position += 1;

    def assign_to_parent(self, chromosome):
        curr_vcf = self.vcf_dfs[chromosome];
        for dnv in self.to_phase:
            self.phased_to_parent[dnv] = [];
            dnv_index = curr_vcf.index[curr_vcf['POS'] == dnv].item();
            de_novo_hap = curr_vcf[self.id][dnv_index];
            for var in self.to_phase[dnv]:
                index = curr_vcf.index[curr_vcf['POS'] == var].item();
                if len(curr_vcf['REF'][index]) > 1 or len(curr_vcf['ALT'][index]) > 1:
                    print("Skipped variant at position " + str(vcf_df['POS'][index]));
                    continue;
                child = curr_vcf[self.id][index];
                mom = curr_vcf[self.id][index];
                dad = curr_vcf[self.id][index];
                if child[:3] == de_novo_hap[:3]:
                    if mom[1:3] == "/1":
                        self.phased_to_parent[dnv].append("mom");
                    else:
                        self.phased_to_parent[dnv].append("dad");
                if child[:3] != de_novo_hap[:3]:
                    if mom[1:3] == "/1":
                        self.phased_to_parent[dnv].append("dad");
                    else:
                        self.phased_to_parent[dnv].append("mom");
