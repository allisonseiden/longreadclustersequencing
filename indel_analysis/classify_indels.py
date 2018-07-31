import pandas as pd
import numpy as np


class Bedfile:

    def __init__(self, bed_location):
        self.orig_bed = pd.read_table(bed_location, names=['Chrom', 'Start', 'End', 'Ref', 'Alt', 'ID'],
                                        sep='\t', engine='python');
        self.mod_bed = self.orig_bed;

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
                if allele_len < 6:
                    self.mod_bed.loc[i, 'Start'] -= 6;
                    self.mod_bed.loc[i, 'End'] += 6;
                else:
                    self.mod_bed.loc[i, 'Start'] -= 2*allele_len;
                    self.mod_bed.loc[i, 'End'] += 2*allele_len;
            else:
                if allele_len < 6:
                    self.mod_bed.loc[i, 'Start'] -= 6;
                    self.mod_bed.loc[i, 'End'] += 6 + allele_len;
                else:
                    self.mod_bed.loc[i, 'Start'] -= 2*allele_len;
                    self.mod_bed.loc[i, 'End'] += 3*allele_len;





if __name__ == '__main__':
    test = Bedfile('/Users/allisonseiden/Documents/longreadclustersequencing/data/1-00801_dnv.bed');
    test.get_indels_from_bed();
    test.get_allele();
    test.change_bounds();
    print(test.mod_bed);
