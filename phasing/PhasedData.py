"""Assign parent of origin.

    ----------------------------------------------------------------------------
    A class to create an object for the phased data for each proband and
    his/her parents

    Each object has variables:
    id - the sample ID for the proband
    mom - the sample ID for the mother
    dad - the sample ID for the father
    bed - a dataframe with info from BED file concerning de novos
    vcf_dfs - a dictionary with key: chromosome number,
                value: dataframe with info from phased vcf file for
                corresponding chromosome
    gtf_dfs - a dictionary with key: a chromosome number,
                value: dataframe with info from gtf file
    dnvs - a dictionary with key: chromosome number, value: list of de novos
            for corresponding chromosome (only ones that whatshap phased)
    bounds - a dictionary with key: chromosome number, value: a dictionary with
                key: de novo position, value: list with discontinuity bounds
    to_phase - a dictionary with key: chromosome number, value: a dictionary
                with key: de novo position, value: list of positions of
                informative variants
    phased_to_parent - a dictionary with key: chromosome number, value: a
                dictionary with key: de novo position, value: list of
                'mom'/'dad' assignments for each informative variant
    parent_df - a dataframe with columns for ID number (ID), chromosome number
                (Chrom), variant position (Location), the number of informative
                reads phased to mom (Mom Count), number of informative reads
                to dad (Dad Count), and columns to mark if the de novo came
                from mom (From Mom), from dad (From Dad), or if it needs a
                second look (Troubleshoot) or is unphased (Unphased)
    ----------------------------------------------------------------------------
"""

import os
import pandas as pd
# import numpy as np


class PhasedData(object):
    def __init__(self, fam_id, trio_df=None, home_dir='/hpc/users/seidea02/'):
        self.id = fam_id
        self.mom = fam_id + '-01'
        self.dad = fam_id + '-02'
        # if WGS VCF ID is different from family ID fill it in here
        if trio_df is not None:
            trio_df_id_row = trio_df.loc[trio_df.Fam_ID == fam_id]
            self.vcf_id = trio_df_id_row['Child'].to_string(index=False)
            self.vcf_id_mom = trio_df_id_row['Mother'].to_string(index=False)
            self.vcf_id_dad = trio_df_id_row['Father'].to_string(index=False)
            # possibly also need to set parent mom and dad IDs
        else:
            self.vcf_id = self.id
            self.vcf_id_mom = self.mom
            self.vcf_id_dad = self.dad
        dnv_f = '{}longreadclustersequencing/data/gmkf2/{}_dnv.bed'.format(
            home_dir, self.id)
        dnv_cols = ['Chrom', 'Start', 'End', 'Ref', 'Var', 'ID']
        self.bed = pd.read_table(dnv_f, sep='\t', names=dnv_cols)
        self.vcf_dfs = {}
        self.gtf_dfs = {}
        self.dnvs = {}
        self.bounds = {}
        self.to_phase = {}
        self.phased_to_parent = {}
        self.parent_df = pd.DataFrame({'ID': [], 'Chrom': [], 'Location': [],
                                       'Mom Count': [], 'Dad Count': [],
                                       'From Mom': [], 'From Dad': [],
                                       'Unphased': []})

    """
        ------------------------------------------------------------------------
        Method to fill vcf_dfs with dataframes created from phased vcfs,
        gtf_dfs with dataframes created from gtf files for indel vcfs, and
        bed file created from first pass of no indels vcfs
        ------------------------------------------------------------------------
    """
    def create_vcf_dictionary(self, whatshap_prefix):
        for i in range(1, 23):
            num = str(i)
            whatshap_vcf = whatshap_prefix.format(
                self.vcf_id, self.id, num) + '.vcf'
            vcf_cols = ['CHROM', 'POS', 'ID', 'REF',
                        'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                        self.vcf_id, self.vcf_id_mom, self.vcf_id_dad]
            if not os.path.exists(whatshap_vcf):
                print('Whatshap VCF not yet made for', whatshap_vcf)
                # Switch to break once program works: only want to run
                # analyses if all data is ready
                continue
            else:
                print('Loading {}...'.format(whatshap_vcf))
            whatshap_gtf = whatshap_prefix.format(
                self.vcf_id, self.id, num) + '.gtf'
            gtf_cols = ['Chrom', 'Allison', 'Start', 'End', 'Felix',
                        'Plus', 'Dot', 'Madeline']
            if not os.path.exists(whatshap_gtf):
                print('Whatshap GTF not yet made for', whatshap_gtf)
                # Switch to break once program works: only want to run
                # analyses if all data is ready
                continue
            self.vcf_dfs["chr{0}".format(i)] = pd.read_table(
                whatshap_vcf, sep='\t', names=vcf_cols, comment='#')
            # for compatibility, replace vcf_ids with ids
            # rename format is old_name:new_name
            self.vcf_dfs["chr{0}".format(i)].rename(
                columns={self.vcf_id: self.id, self.vcf_id_mom: self.mom,
                         self.vcf_id_dad: self.dad},
                inplace=True)
            self.gtf_dfs["chr{0}".format(i)] = pd.read_table(
                whatshap_gtf, sep='\t', names=gtf_cols)
        print('---VCF and GTF dictionaries created for ' + self.id)

    """
        ------------------------------------------------------------------------
        Method to fill vcf_dfs with dataframes created from phased vcfs without
        --indels, gtf_dfs with dataframes created from gtf files for no indel
        vcfs
        ------------------------------------------------------------------------
    """

    def create_vcf_no_indels(self):
        for i in range(1,23):
            num = str(i)
            self.vcf_dfs["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '_no_indels/' + self.id + '_chr' + num + '_phased.vcf',
                                                sep='\t', names=['CHROM', 'POS', 'ID', 'REF',
                                                'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                                self.id, self.mom, self.dad],
                                                comment='#')
            self.gtf_dfs["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '_no_indels/' + self.id + '_chr' + num + '_phased.gtf',
                                                sep='\t', names=['Chrom', 'Allison', 'Start', 'End', 'Felix', 'Plus', 'Dot', 'Madeline'])
        print('---VCF and GTF dictionaries created for ' + self.id)

    """
        ------------------------------------------------------------------------
        Method to fill dnvs dictionary from BED file
        ------------------------------------------------------------------------
    """

    def create_dnvs_dictionary(self):
        num_dnvs = self.bed.shape[0]
        chrom_list = []
        # Ensure that only chromosomes with de novos present in BED file are
        # included and not repeated
        for row in range(0, num_dnvs):
            chr_num = self.bed.loc[row, 'Chrom']
            if chr_num not in chrom_list:
                chrom_list.append(chr_num)

        # Add de novo positions to list, de novo list value of dictionary
        # corresponding to key value of chromosome where de novo is located
        for chrom in chrom_list:
            indices = self.bed.index[self.bed['Chrom'] == chrom].tolist()
            self.dnvs[chrom] = []
            for index in indices:
                self.dnvs[chrom].append(self.bed['End'][index])

        print('---DNV dictionary created for ' + self.id)

    """
        ------------------------------------------------------------------------
        Helper method to search for discontinuities within phased vcf, makes
        use of gtf files to know where upper and lower bounds are
        ------------------------------------------------------------------------
    """

    def search_discon(self, chromosome):
        chr_bounds = {}
        # Collect correct phased VCF file for the current chromosome
        curr_vcf = self.vcf_dfs[chromosome]
        # Collect start and end positions of haplotype blocks from GTF file for
        # current chromosome
        start_list = self.gtf_dfs[chromosome]['Start'].tolist()
        end_list = self.gtf_dfs[chromosome]['End'].tolist()
        # Loop through de novo positions within list of de novos corresponding
        # to current chromosome
        for dnv in self.dnvs[chromosome]:
            chr_bounds[dnv] = []
            # Get index of de novo position within VCF file, will need to look
            # at lines above and below to find discontinuities
            dnv_index = curr_vcf.index[curr_vcf['POS'] == dnv].item()
            # Initially set haplotype to de novo haplotype
            hap = curr_vcf[self.id][dnv_index]
            # If the de novo is unphased, skip over it
            if hap[:3] == "0/1":
                continue
            u_discon = dnv_index
            # Loop through lines around de novo in the VCF file until a variant
            # outside the haplotype block is found, ignore unphased indels
            while (curr_vcf['POS'][u_discon] not in start_list) or (hap[:3] == "0/1" and len(curr_vcf['REF'][u_discon]) > 1) or (hap[:3] == "0/1" and len(curr_vcf['ALT'][u_discon]) > 1):
                u_discon -= 1
                hap = curr_vcf[self.id][u_discon]
            chr_bounds[dnv].append(curr_vcf['POS'][u_discon])
            hap = curr_vcf[self.id][dnv_index]
            l_discon = dnv_index
            distance = abs(dnv - (curr_vcf['POS'][l_discon]))
            while (curr_vcf['POS'][l_discon] not in end_list) or (hap[:3] == "0/1" and len(curr_vcf['REF'][l_discon]) > 1) or (hap[:3] == "0/1" and len(curr_vcf['ALT'][l_discon]) > 1):
                l_discon += 1
                hap = curr_vcf[self.id][l_discon]
            chr_bounds[dnv].append(curr_vcf['POS'][l_discon])
        return chr_bounds

    """
        ------------------------------------------------------------------------
        Method to fill bounds dictionary using search_discon helper method
        ------------------------------------------------------------------------
    """

    def fill_bounds_dictionary(self):
        for chr in self.vcf_dfs:
            self.bounds[chr] = self.search_discon(chr)

        print('---Bounds dictionary created for ' + self.id)

    """
        ------------------------------------------------------------------------
        Helper method to find informative variants between bounds for each dnv
        ------------------------------------------------------------------------
    """

    def find_variants_for_phasing_chr(self, chromosome, n):
        chr_phase = {}
        curr_vcf = self.vcf_dfs[chromosome]
        for dnv in self.dnvs[chromosome]:
            curr_bounds = self.bounds[chromosome][dnv]
            chr_phase[dnv] = []
            # If the bounds list is empty, i.e. the de novo is unphased or
            # there were no variants within the haplotype block, skip
            if len(curr_bounds) == 0:
                continue
            u_index = curr_vcf.index[curr_vcf['POS'] == curr_bounds[0]].item()
            l_index = curr_vcf.index[curr_vcf['POS'] == curr_bounds[1]].item()
            position = u_index
            # Loop through variants between upper and lower bounds and
            # determine which are informative (de novo is phased, and mom
            # and dad do not have the same haplotype)
            while position <= l_index:
                child = curr_vcf[self.id][position]
                mom = curr_vcf[self.mom][position]
                dad = curr_vcf[self.dad][position]
                if child[:3] == "0|1" or child[:3] == "1|0":
                    if mom[:3] != dad[:3]:
                        chr_phase[dnv].append(curr_vcf['POS'][position])
                # Ensure that only the closest n variants are used above and
                # below the de novo
                if (len(chr_phase[dnv]) != 0) and (curr_vcf['POS'][position] < dnv) and (len(chr_phase[dnv]) > n):
                    chr_phase[dnv] = chr_phase[dnv][-n:]
                position += 1
            if len(chr_phase[dnv]) > 2*n:
                chr_phase[dnv] = chr_phase[dnv][:(2*n)]
            else:
                continue
        return chr_phase

    """
        ------------------------------------------------------------------------
        Method to fill to_phase dictionary using find_variants_for_phasing_chr
        helper method
        ------------------------------------------------------------------------
    """

    def find_variants_for_phasing(self, n):
        for chr in self.dnvs:
            self.to_phase[chr] = self.find_variants_for_phasing_chr(chr, n)

        print('---Variants to phase dictionary created for ' + self.id)

    """
        ------------------------------------------------------------------------
        Helper method to conduct the logic for assigning de novos to parent,
        used to keep code clean
        ------------------------------------------------------------------------
    """

    def assign_to_parent_logic(self, index, vcf, dnv_hap, chr_parent, dnv):
        # Collect haplotype information for child, mother, and father from VCF
        child = vcf[self.id][index]
        mom = vcf[self.mom][index]
        dad = vcf[self.dad][index]
        kiddo = child[:3]
        ma = mom[:3]
        pa = dad[:3]
        # Determine which parent the '1' or '0' came from based on whether mom
        # and dad are homozygous or heterozygous for a certain variant and then
        # assign the chromosomes to a parent for each informative variant
        if (ma == '0/1' and pa == '1/1') or (ma == '0/0' and pa == '0/1') or (ma == '0/0' and pa == '1/1'):
            # Compare the informative variants to the de novo (whether both
            # are 0|1 or 1|0 or they are opposite) and assign accordingly
            if kiddo == dnv_hap[:3]:
                chr_parent[dnv].append('dad')
            else:
                chr_parent[dnv].append('mom')
        if ((ma == '0/1') and (pa == '0/0')) or ((ma == '1/1') and (pa == '0/0')) or ((ma == '1/1') and (pa == '0/1')):
            if kiddo == dnv_hap[:3]:
                chr_parent[dnv].append('mom')
            else:
                chr_parent[dnv].append('dad')

    """
        ------------------------------------------------------------------------
        Helper method to assign each informative variant to either mom or dad
        using assign_to_parent_logic helper method
        ------------------------------------------------------------------------
    """

    def assign_to_parent_by_chr(self, chromosome):
        chr_parent = {}
        curr_vcf = self.vcf_dfs[chromosome]
        for dnv in self.to_phase[chromosome]:
            chr_parent[dnv] = []
            dnv_index = curr_vcf.index[curr_vcf['POS'] == dnv].item()
            de_novo_hap = curr_vcf[self.id][dnv_index]
            for var in self.to_phase[chromosome][dnv]:
                index = curr_vcf.index[curr_vcf['POS'] == var].item()
                self.assign_to_parent_logic(index, curr_vcf, de_novo_hap, chr_parent, dnv)
        return chr_parent

    """
        ------------------------------------------------------------------------
        Method to fill phased_to_parent dictionary using assign_to_parent_by_chr
        helper method
        ------------------------------------------------------------------------
    """

    def assign_to_parent(self):
        for chr in self.to_phase:
            self.phased_to_parent[chr] = self.assign_to_parent_by_chr(chr)

        print('---DNVs phased to parent for ' + self.id)


    """
        ------------------------------------------------------------------------
        Method to organize data into final dataframe parent_df
        ------------------------------------------------------------------------
    """

    def convert_to_dataframe(self):
        id_list = []
        chrom_list = []
        location_list = []
        mom_count = []
        dad_count = []
        from_mom = []
        from_dad = []
        unphased = []
        for chr in self.phased_to_parent:
            # For each de novo, go through phased_to_parent list and count
            # how many informative variants were assigned to mom and how many
            # were assigned to dad
            for dnv in self.phased_to_parent[chr]:
                id_list.append(self.id)
                chrom_list.append(chr)
                location_list.append(dnv)
                mom = 0
                dad = 0
                for parent in self.phased_to_parent[chr][dnv]:
                    if parent == 'mom':
                        mom += 1
                    elif parent == 'dad':
                        dad += 1
                    else:
                        continue
                mom_count.append(mom)
                dad_count.append(dad)

        self.parent_df['ID'] = id_list
        self.parent_df['Chrom'] = chrom_list
        self.parent_df['Location'] = location_list
        self.parent_df['Mom Count'] = mom_count
        self.parent_df['Dad Count'] = dad_count

        length = self.parent_df.shape[0]

        for i in range(0, length):
            # Collect the number of informative variants assigned to mom and
            # assigned to dad for each de novo
            ma = self.parent_df['Mom Count'][i]
            pa = self.parent_df['Dad Count'][i]
            # If there are no informative variants for both mom and dad,
            # consider the de novo unphased
            if ma == 0 and pa == 0:
                unphased.append(1)
                from_mom.append(0)
                from_dad.append(0)
            # If 85% or more of the informative variants phased the de novo to
            # mom, consider the de novo to have come from mom
            elif ma/(ma + pa) >= .85:
                from_mom.append(1)
                from_dad.append(0)
                unphased.append(0)
            # If 85% or more of the informative variants phased the de novo to
            # dad, consider the de novo to have come from dad
            elif pa/(ma + pa) >= .85:
                from_mom.append(0)
                from_dad.append(1)
                unphased.append(0)
            # If the informative variants present ambiguous information,
            # consider the de novo unphased
            else:
                from_mom.append(0)
                from_dad.append(0)
                unphased.append(1)

        self.parent_df['From Mom'] = from_mom
        self.parent_df['From Dad'] = from_dad
        self.parent_df['Unphased'] = unphased
        self.parent_df = self.parent_df[['ID', 'Chrom', 'Location', 'From Mom',
                                         'From Dad', 'Unphased']]
        self.parent_df.to_csv(path_or_buf=self.id + '_dataframe.txt', sep='\t',
                              float_format='%g', index=False)

    """
        ------------------------------------------------------------------------
        Calls all the functions that are needed for Pacbio
        data (i.e. no indels)
        ------------------------------------------------------------------------
    """
    def pacbio(self):
        self.create_vcf_no_indels()
        self.create_dnvs_dictionary()
        self.fill_bounds_dictionary()
        self.find_variants_for_phasing(7)
        self.assign_to_parent()
        self.convert_to_dataframe()

    """
        ------------------------------------------------------------------------
        Calls all the functions that are needed for Illumina data
        (i.e. with indels)
        ------------------------------------------------------------------------
    """
    def illumina(self, whatshap_prefix):
        self.create_vcf_dictionary(whatshap_prefix)
        self.create_dnvs_dictionary()
        self.fill_bounds_dictionary()
        self.find_variants_for_phasing(7)
        self.assign_to_parent()
        self.convert_to_dataframe()
