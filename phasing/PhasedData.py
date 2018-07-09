import pandas as pd;
import numpy as np;

"""
    ----------------------------------------------------------------------------
    A class to create an object for the phased data for each proband and his/her
    parents

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
    unphased - dataframe with chromosome number and location of unphased dnvs
    troubleshoot - dataframe with chromosome number and location of troublesome
                    dnvs
    dnvs - a dictionary with key: chromosome number, value: list of de novos for
            corresponding chromosome (only ones that whatshap phased)
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
                reads phased to mom (Mom Count), the number of informative reads
                to dad (Dad Count), and columns to mark whether the de novo came
                from mom (From Mom), from dad (From Dad), or if it needs a
                second look (Troubleshoot) or is unphased (Unphased)
    ----------------------------------------------------------------------------
"""

class PhasedData:
    def __init__(self, patientID):
        self.id = patientID;
        self.mom = patientID + '-01';
        self.dad = patientID + '-02';
        self.bed = pd.DataFrame();
        self.vcf_dfs = {};
        self.dnvs = {};
        self.unphased = pd.DataFrame({'Chrom' : [], 'Location' : []});
        self.trouble = pd.DataFrame({'Chrom' : [], 'Location' : []});
        self.bounds = {};
        self.to_phase = {};
        self.phased_to_parent = {};
        self.parent_df = pd.DataFrame({'ID' : [], 'Chrom' : [], 'Location' : [],
                            'Mom Count' : [], 'Dad Count' : [],
                            'From Mom' : [], 'From Dad' : [], 'Troubleshoot' : [],
                            'Unphased' : []});
        self.gtf_dfs = {};

    """
        ------------------------------------------------------------------------
        Method to fill vcf_dfs with dataframes created from phased vcfs,
        gtf_dfs with dataframes created from gtf files for indel vcfs, and
        bed file created from first pass of no indels vcfs
        ------------------------------------------------------------------------
    """
    def create_vcf_dictionary(self):
        for i in range(1,23):
            num = str(i);
            self.vcf_dfs["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '/' + self.id + '_chr' + num + '_phased.vcf',
                                                sep='\t', names = ['CHROM', 'POS', 'ID', 'REF',
                                                'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                                self.id, self.mom, self.dad],
                                                comment = '#');
            self.gtf_dfs["chr{0}".format(i)] =  pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '/' + self.id + '_chr' + num + 'phased.gtf',
                                                sep='\t', names = ['Chrom', 'Allison', 'Start', 'End', 'Felix', 'Plus', 'Dot', 'Madeline']);
            self.bed = pd.read_table("/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/" + self.id + "/" + self.id + "_no_indels_dnvs.bed",
                                        sep='\t', names = ['Chrom', 'Start', 'End', 'Ref', 'Var', 'ID']);
        print('---VCF and GTF dictionaries created for ' + self.id);

    """
        ------------------------------------------------------------------------
        Method to fill vcf_dfs with dataframes created from phased vcfs without
        --indels, gtf_dfs with dataframes created from gtf files for no indel
        vcfs, and bed file created from original de novo list
        ------------------------------------------------------------------------
    """

    def create_vcf_no_indels(self):
        for i in range(1,23):
            num = str(i);
            self.vcf_dfs["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '_no_indels/' + self.id + '_chr' + num + '_phased.vcf',
                                                sep='\t', names = ['CHROM', 'POS', 'ID', 'REF',
                                                'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                                self.id, self.mom, self.dad],
                                                comment = '#');
            self.gtf_dfs["chr{0}".format(i)] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' + self.id + '_no_indels/' + self.id + '_chr' + num + 'phased.gtf',
                                                sep='\t', names = ['Chrom', 'Allison', 'Start', 'End', 'Felix', 'Plus', 'Dot', 'Madeline']);
            self.bed = pd.read_table('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/'
                                        + self.id + '.hg38.dnv.bed', sep='\t',
                                        names = ['Chrom', 'Start', 'End', 'Ref', 'Var', 'ID']);
        print('---VCF and GTF dictionaries created for ' + self.id);

    """
        ------------------------------------------------------------------------
        Method to fill dnvs dictionary from BED file, ignores de novos that
        were not phased by whatshap and places them into a dataframe
        ------------------------------------------------------------------------
    """

    def create_dnvs_dictionary(self):
        num_dnvs = self.bed.shape[0];
        chrom_list = [];
        for row in range(0, num_dnvs):
            chr_num = self.bed.loc[row, 'Chrom'];
            if chr_num not in chrom_list:
                chrom_list.append(chr_num);

        unphased_chrom = [];
        unphased_loc = [];

        for chrom in chrom_list:
            indices = self.bed.index[self.bed['Chrom'] == chrom].tolist();
            self.dnvs[chrom] = [];
            length = self.vcf_dfs[chrom].shape[0];
            all_unphased = [];
            for i in range(0, length):
                if self.vcf_dfs[chrom][self.id][i][:3] == "0/1":
                    all_unphased.append(self.vcf_dfs[chrom]['POS'][i]);
            for index in indices:
                if self.bed['End'][index] in all_unphased:
                    unphased_chrom.append(chrom);
                    unphased_loc.append(self.bed['End'][index]);
                    # self.unphased.append(self.bed['End'][index]);
                else:
                    self.dnvs[chrom].append(self.bed['End'][index]);

        self.unphased['Chrom'] = unphased_chrom;
        self.unphased['Location'] = unphased_loc;

        print('---DNV dictionary created for ' + self.id);

    """
        ------------------------------------------------------------------------
        Helper method to search for discontinuities within phased vcf, makes
        use of gtf files to know where upper and lower bounds are
        ------------------------------------------------------------------------
    """

    def search_discon(self, chromosome):
        chr_bounds = {};
        curr_vcf = self.vcf_dfs[chromosome];
        start_list = self.gtf_dfs[chromosome]['Start'].tolist();
        end_list = self.gtf_dfs[chromosome]['End'].tolist();
        for dnv in self.dnvs[chromosome]:
            chr_bounds[dnv] = [];
            dnv_index = curr_vcf.index[curr_vcf['POS'] == dnv].item();
            hap = curr_vcf[self.id][dnv_index];
            u_discon = dnv_index;
            # distance = abs(dnv - (curr_vcf['POS'][u_discon]));
            # while (hap[:3] != "0/1" or (hap[:3] == "0/1" and len(curr_vcf['REF'][u_discon]) > 1) or (hap[:3] == "0/1" and len(curr_vcf['ALT'][u_discon]) > 1)) and distance <= 10000:
            while (curr_vcf['POS'][u_discon] not in start_list) or (hap[:3] == "0/1" and len(curr_vcf['REF'][u_discon]) > 1) or (hap[:3] == "0/1" and len(curr_vcf['ALT'][u_discon]) > 1):
                u_discon -= 1;
                hap = curr_vcf[self.id][u_discon];
                # distance = abs(dnv - (curr_vcf['POS'][u_discon]))
            chr_bounds[dnv].append(curr_vcf['POS'][u_discon]);
            hap = curr_vcf[self.id][dnv_index];
            l_discon = dnv_index;
            distance = abs(dnv - (curr_vcf['POS'][l_discon]));
            # while (hap[:3] != "0/1" or (hap[:3] == "0/1" and len(curr_vcf['REF'][l_discon]) > 1) or (hap[:3] == "0/1" and len(curr_vcf['ALT'][l_discon]) > 1)) and distance <= 10000:
            while (curr_vcf['POS'][l_discon] not in end_list) or (hap[:3] == "0/1" and len(curr_vcf['REF'][l_discon]) > 1) or (hap[:3] == "0/1" and len(curr_vcf['ALT'][l_discon]) > 1):
                l_discon += 1;
                hap = curr_vcf[self.id][l_discon];
                # distance = abs(dnv - (curr_vcf['POS'][l_discon]));
            chr_bounds[dnv].append(curr_vcf['POS'][l_discon]);
        return chr_bounds;

    """
        ------------------------------------------------------------------------
        Method to fill bounds dictionary using search_discon helper method
        ------------------------------------------------------------------------
    """

    def fill_bounds_dictionary(self):
        for chr in self.dnvs:
            self.bounds[chr] = self.search_discon(chr);

        print('---Bounds dictionary created for ' + self.id);

    """
        ------------------------------------------------------------------------
        Helper method to find informative variants between bounds for each dnv
        ------------------------------------------------------------------------
    """

    def find_variants_for_phasing_chr(self, chromosome):
        chr_phase = {};
        curr_vcf = self.vcf_dfs[chromosome];
        for dnv in self.dnvs[chromosome]:
            curr_bounds = self.bounds[chromosome][dnv];
            chr_phase[dnv] = [];
            u_index = curr_vcf.index[curr_vcf['POS'] == curr_bounds[0]].item();
            l_index = curr_vcf.index[curr_vcf['POS'] == curr_bounds[1]].item();
            position = u_index;
            while position <= l_index:
                child = curr_vcf[self.id][position];
                mom = curr_vcf[self.mom][position];
                dad = curr_vcf[self.dad][position];
                if child[:3] == "0|1" or child[:3] == "1|0":
                    if mom[:3] != dad[:3]:
                        chr_phase[dnv].append(curr_vcf['POS'][position]);
                # if (len(chr_phase[dnv]) != 0) and (curr_vcf['POS'][position] < dnv) and (len(chr_phase[dnv]) > n):
                #     chr_phase[dnv] = chr_phase[dnv][-n:]
                position += 1;
            # if len(chr_phase[dnv]) > 2*n:
            #     chr_phase[dnv] = chr_phase[dnv][:(2*n)];
            # else:
            #     continue;
        return chr_phase;

    """
        ------------------------------------------------------------------------
        Method to fill to_phase dictionary using find_variants_for_phasing_chr
        helper method
        ------------------------------------------------------------------------
    """

    def find_variants_for_phasing(self):
        for chr in self.dnvs:
            self.to_phase[chr] = self.find_variants_for_phasing_chr(chr);

        print('---Variants to phase dictionary created for ' + self.id);

    """
        ------------------------------------------------------------------------
        Helper method to conduct the logic for assigning de novos to parent,
        used to keep code clean
        ------------------------------------------------------------------------
    """

    def assign_to_parent_logic(self, index, vcf, dnv_hap, chr_parent, dnv):
        child = vcf[self.id][index];
        mom = vcf[self.mom][index];
        dad = vcf[self.dad][index];
        kiddo = child[:3];
        ma = mom[:3];
        pa = dad[:3];
        if (ma == '0/1' and pa == '1/1') or (ma == '0/0' and pa == '0/1') or (ma == '0/0' and pa == '1/1'):
            if kiddo == dnv_hap[:3]:
                chr_parent[dnv].append('dad');
            else:
                chr_parent[dnv].append('mom');
        if ((ma == '0/1') and (pa == '0/0')) or ((ma == '1/1') and (pa == '0/0')) or ((ma == '1/1') and (pa == '0/1')):
            if kiddo == dnv_hap[:3]:
                chr_parent[dnv].append('mom');
            else:
                chr_parent[dnv].append('dad');

    """
        ------------------------------------------------------------------------
        Helper method to assign each informative variant to either mom or dad
        using assign_to_parent_logic helper method
        ------------------------------------------------------------------------
    """

    def assign_to_parent_by_chr(self, chromosome):
        chr_parent = {};
        curr_vcf = self.vcf_dfs[chromosome];
        for dnv in self.to_phase[chromosome]:
            chr_parent[dnv] = [];
            dnv_index = curr_vcf.index[curr_vcf['POS'] == dnv].item();
            de_novo_hap = curr_vcf[self.id][dnv_index];
            for var in self.to_phase[chromosome][dnv]:
                index = curr_vcf.index[curr_vcf['POS'] == var].item();
                self.assign_to_parent_logic(index, curr_vcf, de_novo_hap, chr_parent, dnv);
        return chr_parent;

    """
        ------------------------------------------------------------------------
        Method to fill phased_to_parent dictionary using assign_to_parent_by_chr
        helper method
        ------------------------------------------------------------------------
    """

    def assign_to_parent(self):
        for chr in self.to_phase:
            self.phased_to_parent[chr] = self.assign_to_parent_by_chr(chr);

        print('---DNVs phased to parent for ' + self.id);


    """
        ------------------------------------------------------------------------
        Method to organize data into final dataframe parent_df
        ------------------------------------------------------------------------
    """

    def convert_to_dataframe(self):
        id_list = [];
        chrom_list = [];
        location_list = [];
        mom_count = [];
        dad_count = [];
        from_mom = [];
        from_dad = [];
        trouble = [];
        unphased = [];
        for chr in self.phased_to_parent:
            for dnv in self.phased_to_parent[chr]:
                id_list.append(self.id);
                chrom_list.append(chr);
                location_list.append(dnv);
                mom = 0;
                dad = 0;
                for parent in self.phased_to_parent[chr][dnv]:
                    if parent == 'mom':
                        mom += 1;
                    elif parent == 'dad':
                        dad += 1;
                    else:
                        continue;
                mom_count.append(mom);
                dad_count.append(dad);

        self.parent_df['ID'] = id_list;
        self.parent_df['Chrom'] = chrom_list;
        self.parent_df['Location'] = location_list;
        self.parent_df['Mom Count'] = mom_count;
        self.parent_df['Dad Count'] = dad_count;

        length = self.parent_df.shape[0];

        for i in range(0, length):
            ma = self.parent_df['Mom Count'][i];
            pa = self.parent_df['Dad Count'][i];
            if ma == 0 and pa == 0:
                trouble.append(1);
                from_mom.append(0);
                from_dad.append(0);
            elif ma/(ma + pa) >= .85:
                from_mom.append(1);
                from_dad.append(0);
                trouble.append(0);
            elif pa/(ma + pa) >= .85:
                from_mom.append(0);
                from_dad.append(1);
                trouble.append(0);
            else:
                from_mom.append(0);
                from_dad.append(0);
                trouble.append(1);

        self.parent_df['From Mom'] = from_mom;
        self.parent_df['From Dad'] = from_dad;
        self.parent_df['Troubleshoot'] = trouble;


        self.parent_df = self.parent_df[['ID', 'Chrom', 'Location', 'Mom Count',
                                            'Dad Count', 'From Mom', 'From Dad',
                                            'Troubleshoot']];

        trouble_chrom = [];
        trouble_loc = [];
        for i in range(0, length):
            if self.parent_df['Troubleshoot'][i] == 1:
                trouble_chrom.append(self.parent_df['Chrom'][i]);
                trouble_loc.append(self.parent_df['Location'][i]);

        self.trouble['Chrom'] = trouble_chrom;
        self.trouble['Location'] = trouble_loc;

        self.parent_df = self.parent_df.groupby('ID').sum();
        self.parent_df['Unphased'] = self.unphased.shape[0];
        self.parent_df = self.parent_df.loc[:,['From Mom', 'From Dad', 'Troubleshoot', 'Unphased']];

    """
        ------------------------------------------------------------------------
        Writes troubleshoot and unphased de novo locations to bed file to pass
        to indels vcfs
        ------------------------------------------------------------------------
    """


    def write_to_bed(self):
        filename = "/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/" + self.id + "/" + self.id + "_no_indels_dnvs.bed";
        bed_file = open(filename, "w");
        trouble_length = self.trouble.shape[0];
        for i in range(0, trouble_length):
            line = str(self.trouble['Chrom'][i]) + '\t.\t' + str(self.trouble['Location'][i]) + '\t.\t.\t' + self.id + '\n';
            bed_file.write(line);
        unphased_length = self.unphased.shape[0];
        for i in range(0, unphased_length):
            line = str(self.unphased['Chrom'][i]) + '\t.\t' + str(self.unphased['Location'][i]) + '\t.\t.\t' + self.id + '\n';
            bed_file.write(line);
        bed_file.close();
