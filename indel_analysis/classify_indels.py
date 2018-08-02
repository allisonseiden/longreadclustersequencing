import pandas as pd
import numpy as np
import subprocess as sp
import argparse

"""
    Sorts mutations that are insertions or deletions into mutational classes:
    HR: homopolymer run (mutation is in a region where there are 6 or more
    copies of the nucleotide being inserted or deleted)
    CCC: change in copy count (the allele being inserted or deleted has 1 or
    more repeats in the mutation region)
    non-CCC: no change in copy count (the allele being inserted or deleted is
    not repeated in the mutation region)

    Arguments:
    --bed       BED file with all variants. Can include multiple IDs or a
                single ID but must be formatted as Chrom-Start-End-Ref-Alt-ID
    --fasta     reference fasta file
    --repeat    RepeatMasker file downloaded from UCSC Genome Browser,
                instructions on how to download can be found in the docs

    Written by: Allison Seiden
                ahseiden@gmail.com

"""
class SortIt:

    def __init__(self, bed_location, fasta_location, repeat_masker):
        self.orig_bed = pd.read_table(bed_location, names=['Chrom', 'Start',
                                        'End', 'Ref', 'Alt', 'ID'],
                                        sep='\t', engine='python');
        self.mod_bed = self.orig_bed;
        self.fasta = fasta_location;
        self.repeat_masker = repeat_masker;

    """ Retrieves indels from BED file of all variants
    """
    def get_indels_from_bed(self):
        length = self.mod_bed.shape[0];
        indices = [];
        for i in range(length):
            ref_len = len(self.mod_bed.loc[i, 'Ref']);
            alt_len = len(self.mod_bed.loc[i, 'Alt']);
            if ref_len == 1 and alt_len == 1:
                indices.append(i);
        self.mod_bed = self.mod_bed.drop(self.mod_bed.index[indices]);
        self.mod_bed.reset_index(inplace=True);
        self.mod_bed = self.mod_bed[['Chrom', 'Start', 'End', 'Ref', 'Alt', 'ID']];

    """ Collects base sequence that is inserted or deleted from Ref/Alt
    """
    def get_allele(self):
        length = self.mod_bed.shape[0];
        for i in range(length):
            if len(self.mod_bed.loc[i, 'Ref']) < len(self.mod_bed.loc[i, 'Alt']):
                self.mod_bed.loc[i, 'Allele'] = self.mod_bed.loc[i, 'Alt'][1:];
            else:
                self.mod_bed.loc[i, 'Allele'] = self.mod_bed.loc[i, 'Ref'][1:];

    """ Changes Start/End locations in modified BED file in order to use
        bedtools getfasta to get surrounding sequence
    """
    def change_bounds(self):
        length = self.mod_bed.shape[0];
        for i in range(length):
            allele_len = len(self.mod_bed.loc[i, 'Allele']);
            # Insertions
            if len(self.mod_bed.loc[i, 'Ref']) < len(self.mod_bed.loc[i, 'Alt']):
                if allele_len == 1:
                    self.mod_bed.loc[i, 'Start'] -= 6;
                    self.mod_bed.loc[i, 'End'] += 6;
                else:
                    self.mod_bed.loc[i, 'Start'] -= 2*allele_len;
                    self.mod_bed.loc[i, 'End'] += 2*allele_len;
            # Deletions
            else:
                if allele_len == 1:
                    self.mod_bed.loc[i, 'Start'] -= 5;
                    self.mod_bed.loc[i, 'End'] += 7;
                else:
                    self.mod_bed.loc[i, 'Start'] -= (2*allele_len-1);
                    self.mod_bed.loc[i, 'End'] += 3*allele_len;

    """ Uses bedtools getfasta to get sequence surrounding indel
    """
    def get_fasta(self):
        self.mod_bed.to_csv(path_or_buf='tmp.bed', sep='\t', header=False,
                            index=False);
        cmd = 'bedtools getfasta -fi ' + self.fasta + ' -bed tmp.bed -fo ';
        cmd += 'fasta_tmp.bed -tab';
        sp.call(cmd, shell=True);
        fasta_df = pd.read_table('fasta_tmp.bed', sep=':|-|\t', engine='python',
                                    names=['Chrom', 'Start', 'End', 'Sequence']);
        self.mod_bed.set_index(['Chrom', 'Start', 'End'], inplace=True);
        fasta_df.set_index(['Chrom', 'Start', 'End'], inplace=True);
        self.mod_bed = self.mod_bed.join(fasta_df, how='left');
        self.mod_bed.reset_index(inplace=True);
        self.mod_bed['Sequence'] = self.mod_bed['Sequence'].str.upper();
        sp.call('rm tmp.bed fasta_tmp.bed', shell=True);


    """ Retrieves the bases adjacent to indel.
        Number of bases retrieved on either side of indel depends on allele
        length
        Returns tuple with previous bases and next bases as elements
    """
    def get_prev_next_bases(self):
        length = self.mod_bed.shape[0];
        prev_next_bases = [];
        for i in range(length):
            allele = self.mod_bed.loc[i, 'Allele'];
            allele_len = len(self.mod_bed.loc[i, 'Allele']);
            seq = self.mod_bed.loc[i, 'Sequence'];
            mid = int(len(seq)/2);
            half_allele = int(allele_len/2);
            seq_bef = "";
            seq_aft = "";
            ind = mid + 1;
            # Insertions
            if len(self.mod_bed.loc[i, 'Ref']) < len(self.mod_bed.loc[i, 'Alt']):
                seq_bef = seq[mid-allele_len+1:mid+1];
                # handles indexing error of even-length and odd-length alleles
                # behaving differently
                if len(seq) % 2 != 0:
                    seq_aft = seq[mid+1:mid+allele_len+1];
                else:
                    seq_aft = seq[mid:mid+allele_len];
            # Deletions
            else:
                seq_bef = seq[(mid-half_allele-allele_len):(mid-half_allele)];
                if len(seq) % 2 != 0:
                    seq_aft = seq[(mid+half_allele+1):(mid+half_allele+allele_len+1)];
                else:
                    seq_aft = seq[(mid+half_allele):(mid+half_allele+allele_len)];
            # check for homopolymer or change in copy count of single bases
            if allele_len == 1:
                seq_aft = seq[ind:ind+1];
                while (ind < len(seq)) and (seq[ind:ind+1] == allele):
                    seq_aft += seq[ind:ind+1];
                    ind += 1;
            prev_next_bases.append([seq_bef, seq_aft]);
        self.mod_bed = self.mod_bed[['Chrom', 'Start', 'End', 'Ref',
                                        'Alt', 'Allele', 'ID']]
        return prev_next_bases;

    """ Assign mutational class using sequences adjacent to indel
    """
    def assign_class(self):
        prev_next_bases = self.get_prev_next_bases();
        length = self.mod_bed.shape[0];
        for i in range(length):
            seq_before = prev_next_bases[i][0];
            seq_after = prev_next_bases[i][1];
            allele = self.mod_bed.loc[i, 'Allele'];
            if len(allele) == 1:
                # homopolymer
                if len(seq_before) >= 6 or len(seq_after) >= 6:
                    self.mod_bed.loc[i, 'Indel_Class'] = 'HR';
                # change in copy count
                elif (len(seq_before) > 1 or len(seq_after) > 1 or
                    seq_before == allele or seq_after == allele):
                    self.mod_bed.loc[i, 'Indel_Class'] = 'CCC';
                # no change in copy count
                else:
                    self.mod_bed.loc[i, 'Indel_Class'] = 'non-CCC';
            else:
                # change in copy count
                if allele == seq_before or allele == seq_after:
                    self.mod_bed.loc[i, 'Indel_Class'] = 'CCC';
                # no change in copy count
                else:
                    self.mod_bed.loc[i, 'Indel_Class'] = 'non-CCC';

    def intersect_repeat(self):
        self.mod_bed.to_csv(path_or_buf='tmp.bed', sep='\t', header=False,
                            index=False);
        cmd = 'bedtools intersect -a tmp.bed -b ' + self.repeat_masker;
        cmd += ' -wb -loj > tmp_intersect.bed';
        sp.call(cmd, shell=True);
        repeat_df = pd.read_table('tmp_intersect.bed', sep='\t',
                                    names=['Chrom', 'Start', 'End', 'Ref',
                                    'Alt', 'Allele', 'ID', 'Indel_Class',
                                    'genoName', 'genoStart', 'genoEnd',
                                    'repName', 'repClass', 'repFamily']);
        sp.call('rm tmp.bed tmp_intersect.bed', shell=True);

        self.mod_bed.set_index(['Chrom', 'Start', 'End', 'Ref', 'Alt',
                                'Allele', 'ID', 'Indel_Class'], inplace=True);
        repeat_df.set_index(['Chrom', 'Start', 'End', 'Ref', 'Alt',
                                'Allele', 'ID', 'Indel_Class'], inplace=True);
        self.mod_bed = self.mod_bed.join(repeat_df, how='left');
        self.mod_bed.reset_index(inplace=True);
        self.mod_bed = self.mod_bed[['ID', 'Chrom', 'Ref',
                                        'Alt', 'Allele', 'Indel_Class', 'repName',
                                        'repClass', 'repFamily']];

        # reassign start and end columns to original locations
        self.orig_bed.set_index(['ID', 'Chrom', 'Ref', 'Alt'], inplace=True);
        # temp = self.orig_bed.join(self.mod_bed, how='left');
        print(self.orig_bed);
        # self.mod_bed['Start'] = self.orig_bed['Start'];
        # self.mod_bed['End'] = self.orig_bed['End'];

def main():
    parser = argparse.ArgumentParser(description="Sorts indels into " +
                                        "mutational classes", epilog="Allison" +
                                        " Seiden <ahseiden@gmail.com>");
    parser.add_argument("-b", "--bed", help="Location of BED file with all " +
                        "variants. Must be formatted as " +
                        "Chrom/Start/End/Ref/Alt/PatientID.", required=True);
    parser.add_argument("-f", "--fasta", help="Location of reference fasta file.",
                        required=True);
    parser.add_argument("-r", "--repeat", help="Location of RepeatMasker file " +
                        "downloaded from UCSC Genome Browser. Refer to docs " +
                        "to see how to download RepeatMasker.", required=True);
    args = parser.parse_args();

    ravenclaw = SortIt(args.bed, args.fasta, args.repeat);
    ravenclaw.get_indels_from_bed();
    ravenclaw.get_allele();
    ravenclaw.change_bounds();
    ravenclaw.get_fasta();
    ravenclaw.assign_class();
    ravenclaw.intersect_repeat();
    # ravenclaw.mod_bed.to_csv(path_or_buf='classified_indels.txt', sep='\t', header=False, index=False);

if __name__ == '__main__':
    main();
