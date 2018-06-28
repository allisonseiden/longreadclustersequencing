import pandas as pd
import numpy as np

chr22_100801 = pd.read_table('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/1-00801/1-00801_chr22_phased.vcf',
                                sep='\t', names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '1-00801', '1-00801-01', '1-00801-02'],
                                comment = '#');

#print(chr22_100801.loc[:,['POS','1-00801']]);
pos_hap = pd.DataFrame(chr22_100801.loc[:['POS', '1-00801']]);
test = pd.Series(pos_hap.loc[:,'1-00801']);
print(test.size);
