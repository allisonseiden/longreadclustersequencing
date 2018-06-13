import pandas as pd
import numpy as np

patient_100801 = pd.read_table('~/Documents/Summer2018/bed_files/1-00801.hg38.dnv.bed',
                                sep='\t', names = ['Chr', 'Start', 'End', 'Ref', 'Var', 'ID'])
#print(patient_100801)
grouped_by_chr = patient_100801.groupby('Chr')
