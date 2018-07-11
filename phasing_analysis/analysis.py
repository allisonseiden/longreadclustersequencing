import pandas as pd
import numpy as np

# last patient ID not in list because fixing bam
patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389', '1-04460',
                '1-04537', '1-05443', '1-05673'];
bed_dictionary = {};
dnv_df = pd.DataFrame();

for ID in patientIDs:
    bed_dictionary[ID] = pd.read_table('/hpc/users/seidea02/www/PacbioProject/DNV_calls/BED/' + ID + '.hg38.dnv.bed',
                                sep='\t', names = ['Chrom', 'Start', 'Loc', 'Ref', 'Alt', 'ID']);
    dnv_df.append(bed_dictionary[ID]);

print(dnv_df);
