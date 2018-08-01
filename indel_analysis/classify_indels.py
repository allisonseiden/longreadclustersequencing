import pandas as pd
import numpy as np
import subprocess as sp


class Bedfile:

    def __init__(self, bed_location, fasta_location):
        self.orig_bed = pd.read_table(bed_location, names=['Chrom', 'Start', 'End', 'Ref', 'Alt', 'ID'],
                                        sep='\t', engine='python');
        self.mod_bed = self.orig_bed;
        self.fasta = fasta_location;

    def get_indels_from_bed(self):
        length = self.mod_bed.shape[0];
        indices = [];
        for i in range(length):
            if len(self.mod_bed.loc[i, 'Ref']) == 1 and len(self.mod_bed.loc[i, 'Alt']) == 1:
                indices.append(i);
        self.mod_bed = self.mod_bed.drop(self.mod_bed.index[indices]);
        self.mod_bed.reset_index(inplace=True);
        self.mod_bed = self.mod_bed[['Chrom', 'Start', 'End', 'Ref', 'Alt']];

    def get_allele(self):
        # put in error handling for if their ref alt is all messed up
        length = self.mod_bed.shape[0];
        for i in range(length):
            if len(self.mod_bed.loc[i, 'Ref']) == 1:
                self.mod_bed.loc[i, 'Allele'] = self.mod_bed.loc[i, 'Alt'][1:];
            else:
                self.mod_bed.loc[i, 'Allele'] = self.mod_bed.loc[i, 'Ref'][1:];

    def change_bounds(self):
        length = self.mod_bed.shape[0];
        for i in range(length):
            allele_len = len(self.mod_bed.loc[i, 'Allele']);
            if len(self.mod_bed.loc[i, 'Ref']) < len(self.mod_bed.loc[i, 'Alt']):
                if allele_len == 1:
                    self.mod_bed.loc[i, 'Start'] -= 6;
                    self.mod_bed.loc[i, 'End'] += 6;
                else:
                    self.mod_bed.loc[i, 'Start'] -= 2*allele_len;
                    self.mod_bed.loc[i, 'End'] += 2*allele_len;
            else:
                if allele_len == 1:
                    self.mod_bed.loc[i, 'Start'] -= 5;
                    self.mod_bed.loc[i, 'End'] += 7;
                else:
                    self.mod_bed.loc[i, 'Start'] -= (2*allele_len-1);
                    self.mod_bed.loc[i, 'End'] += 3*allele_len;

    def get_fasta(self):
        self.mod_bed.to_csv(path_or_buf='tmp.bed', sep='\t', header=False, index=False);
        cmd = 'bedtools getfasta -fi ' + self.fasta + ' -bed tmp.bed -fo fasta_tmp.bed -tab';
        sp.call(cmd, shell=True);
        fasta_df = pd.read_table('fasta_tmp.bed', sep=':|-|\t', engine='python', names=['Chrom', 'Start', 'End', 'Sequence']);
        self.mod_bed.set_index(['Chrom', 'Start', 'End'], inplace=True);
        fasta_df.set_index(['Chrom', 'Start', 'End'], inplace=True);
        self.mod_bed = self.mod_bed.join(fasta_df, how='left');
        self.mod_bed.reset_index(inplace=True);
        sp.call('rm tmp.bed fasta_tmp.bed', shell=True);

    def get_prev_next_bases(self):
        length = self.mod_bed.shape[0];
        for i in range(length):
            allele_len = len(self.mod_bed.loc[i, 'Allele']);
            ref_len = len(self.mod_bed.loc[i, 'Ref']);
            alt_len = len(self.mod_bed.loc[i, 'Alt']);
            seq = self.mod_bed.loc[i, 'Sequence'];
            mid = int(len(seq)/2);
            half_allele = int(allele_len/2);
            bases_before = "";
            bases_after = "";
            index = mid + 1;
            if ref_len < alt_len: # insertions
                bases_before = seq[mid-allele_len+1:mid+1]; # start index at place where insertion occurs (A -> AT includes A in before bases)
                if allele_len == 1: # check for homopolymer or change in copy count of single bases
                    bases_after += seq[index:index+1]; # makes sure that base after is included even if it is not a copy of the allele
                    index += 1;
                    while (index < len(seq)) and (seq[index:index+1] == self.mod_bed.loc[i, 'Allele']):
                        bases_after += seq[index:index+1];
                        index += 1;
                elif len(seq) % 2 != 0: # handles indexing error of even-length and odd-length alleles behaving differently
                    bases_after = seq[mid+1:mid+allele_len+1];
                else:
                    bases_after = seq[mid:mid+allele_len];
            else: # deletions
                bases_before = seq[(mid-half_allele-allele_len):(mid-half_allele)];
                if allele_len == 1:
                    bases_after += seq[index:index+1]; # makes sure that base after is included even if it is not a copy of the allele
                    index += 1;
                    while (index < len(seq)) and (seq[index:index+1] == self.mod_bed.loc[i, 'Allele']):
                        bases_after += seq[index:index+1];
                        index += 1;
                elif len(seq) % 2 != 0:
                    bases_after = seq[(mid+half_allele+1):(mid+half_allele+allele_len+1)];
                else:
                    bases_after = seq[(mid+half_allele):(mid+half_allele+allele_len)];
            print(self.mod_bed.loc[i, 'Allele']);
            print(seq);
            print(bases_before);
            print(bases_after);





if __name__ == '__main__':
    test = Bedfile('/hpc/users/seidea02/longreadclustersequencing/data/1-03897_dnv.bed', '/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa');
    test.get_indels_from_bed();
    test.get_allele();
    test.change_bounds();
    test.get_fasta();
    test.get_prev_next_bases();
